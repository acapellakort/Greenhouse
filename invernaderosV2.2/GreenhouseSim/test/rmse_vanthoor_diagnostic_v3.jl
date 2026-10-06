# =============================================================================
# rmse_vanthoor_diagnostic.jl  (v3 — Tpipe heating estimate + CO2 rate fix)
#
# Changes vs. v2:
#   [D] Tpipe nighttime estimate: when PipeGrow/PipeLow < 5°C (heating off
#       in logger), instead of Tout+2, use:
#         if Tout < T_heat_setpoint - 2°C  →  Tpipe = T_heat_setpoint + PIPE_DELTA_K
#         else                              →  Tpipe = Tout + 2°C
#       where T_heat_setpoint comes from zone-averaged HeatTemp_Vip and
#       PIPE_DELTA_K = 15 K (pipe runs ~15 K above air setpoint when active).
#       Rationale: when outdoor is cold enough that heating is needed, the pipe
#       was almost certainly warm even if the logger shows 0.
#
#   [E] CO2 injection rate: replace bang-bang U10=1 with U10=CO2_INJECT_FRAC
#       (default 0.05) when dosing is active. U10=1 at full model capacity
#       overestimates the real injection rate by ~20×, causing 0–3000 ppm
#       oscillations. CO2_INJECT_FRAC is a tunable parameter.
#
# All v2 fixes remain:
#   [A] CO2 dosing gated by CO2_Vip vs 440 ppm (PPM units corrected)
#   [B] Pipe clamp: readings ≤ 5°C → fallback (now [D] above, not Tout+2)
#   [C] Per-day breakdown table printed
#
# Run from the GreenhouseSim folder:
#   julia --project=. test/rmse_vanthoor_diagnostic_v3.jl
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Statistics, Printf, Plots

# ── paths ────────────────────────────────────────────────────────────────────
const INV     = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TEAM    = "AiCU"
const GH_CSV  = joinpath(INV, "invernaderos2.0", "data", TEAM, "Greenhouse_climate.csv")
const MET_CSV = joinpath(INV, "invernaderos2.0", "data", "meteo.csv")
const VIP_CSV = joinpath(INV, "invernaderos2.0", "data", TEAM, "vip.csv")
const CLM_JS  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
const CRP_JS  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
const FIGDIR  = joinpath(@__DIR__, "figures"); mkpath(FIGDIR)

# ── settings ─────────────────────────────────────────────────────────────────
const DAY0_SERIAL = 43326.0
const N_DAYS      = 42
const LAI_FIX     = 2.0
const EXCEL_UNIX0 = 25569.0
const SIGMA_SB    = 5.670374419e-8

# [A] CO2 dosing gate (PPM)
const CO2_AMB_PPM    = 390.0   # outdoor ambient CO₂ [ppm]
const CO2_THR_PPM    = 50.0    # dose when setpoint > ambient + this [ppm]

# [D] Tpipe heating estimate
const PIPE_DELTA_K   = 25.0    # pipe above air setpoint when heating circuit runs [K]
const PIPE_MIN_VALID = 5.0     # readings ≤ this are treated as "off" [°C]
const PIPE_WARM_TOUT = 22.0   # °C — if T_out exceeds this, never assume heating is running
                              # (real growers don't heat when it's warm outside)

# ── calibrated parameters (t-walk MCMC, Note 28) ─────────────────────────────
# Set USE_CALIBRATED = false to run with engine defaults (sensitivity check).
const USE_CALIBRATED  = true
const CAL_NU4         = 4.46e-05   # leakage [m s⁻¹]   (twalk28, pre-cover-node)
const CAL_GAMMA4      = 44.0       # stomatal resist.   (twalk28, pre-cover-node)
const CAL_K_GROUND    = 24.5       # floor coupling [W/(m²·K)] (twalk28, pre-cover-node; compensatory — re-calibrate after fixing spike days)
const CAL_QOI         = ["nu4", "gamma4", "k_ground"]
const CAL_X           = USE_CALIBRATED ? [CAL_NU4, CAL_GAMMA4, CAL_K_GROUND] : Float64[]
const CAL_NAMES       = USE_CALIBRATED ? CAL_QOI : String[]

# [E] CO2 injection fraction (tune this: 1.0 = full model capacity → 20× too high)
const CO2_INJECT_FRAC = 0.3   # fraction of model's U10 capacity to use

serial2unix(s) = (s - EXCEL_UNIX0) * 86400.0
_num(x) = (x isa Number && !(x isa Missing) && !isnan(Float64(x))) ? Float64(x) : NaN
stepfun(ts, vs) = t -> @inbounds vs[clamp(searchsortedlast(ts, t), 1, length(vs))]
_q2(T) = T > 0 ?
    611.21 * exp((18.678 - T/234.5) * (T / (257.14 + T))) :
    611.21 * exp((23.036 - T/333.7) * (T / (279.82 + T)))
_vpout(Tc, RH) = _q2(Tc) * clamp(RH, 1.0, 100.0) / 100.0

# ── load data ────────────────────────────────────────────────────────────────
gh  = DataFrame(CSV.File(GH_CSV))
mt  = DataFrame(CSV.File(MET_CSV))
vip = DataFrame(CSV.File(VIP_CSV))
n   = min(nrow(gh), nrow(mt), nrow(vip))

win_lo = DAY0_SERIAL
win_hi = DAY0_SERIAL + N_DAYS

ts       = Float64[]
U1v      = Float64[];  U8v     = Float64[];  U10v    = Float64[]
U11v     = Float64[];  U12v    = Float64[];  Tpipev  = Float64[]
Ig_v     = Float64[];  Tout_v  = Float64[];  Wind_v  = Float64[]
VPout_v  = Float64[];  Tsky_v  = Float64[]
Tair_obs = Float64[];  CO2_obs = Float64[];  RH_obs  = Float64[]

# counters for diagnostics
n_pipe_heated = 0   # rows where fallback used heating-setpoint estimate
n_pipe_cold   = 0   # rows where fallback used Tout+2 (warm outdoor, no heating)
n_pipe_meas   = 0   # rows where actual measurement used

for i in 1:n
    g = _num(gh.GHtime[i])
    (isnan(g) || g < win_lo || g >= win_hi) && continue

    # ── outdoor weather ──────────────────────────────────────────────────────
    To  = _num(mt.Tout[i]);   Ig  = _num(mt.Iglob[i])
    Ws  = _num(mt.Windsp[i]); Rho = _num(mt.Rhout[i]);  Py = _num(mt.Pyrgeo[i])
    ToK  = (isnan(To) ? 15.0 : To) + 273.15
    Igv  = isnan(Ig) ? 0.0 : max(Ig, 0.0)
    Wsv  = isnan(Ws) ? 2.0 : Ws
    RHov = isnan(Rho) ? 80.0 : Rho
    Ldn  = SIGMA_SB * ToK^4 + (isnan(Py) ? 0.0 : Py)
    Tskv = isnan(Py) ? (-0.4 + 273.15) : clamp((Ldn / SIGMA_SB)^0.25, ToK-30.0, ToK+2.0)

    # ── actuator states ──────────────────────────────────────────────────────
    vl = _num(gh.VentLee[i]);   vw  = _num(gh.Ventwind[i])
    en = _num(gh.EnScr[i]);     al  = _num(gh.AssimLight[i])
    pg = _num(gh.PipeGrow[i]);  pl  = _num(gh.PipeLow[i])

    U8  = clamp(((isnan(vl) ? 0.0 : vl) + (isnan(vw) ? 0.0 : vw)) / 2.0 / 100.0, 0.0, 1.0)
    U1  = isnan(en) ? 0.0 : clamp(en  / 100.0, 0.0, 1.0)
    U12 = isnan(al) ? 0.0 : clamp(al  / 100.0, 0.0, 1.0)

    # U11: heating setpoint from vip.csv (zone-averaged)
    h1 = _num(vip.HeatTemp_Vip_1[i]); h2 = _num(vip.HeatTemp_Vip_2[i])
    h3 = _num(vip.HeatTemp_Vip_3[i]); h5 = _num(vip.HeatTemp_Vip_5[i])
    h6 = _num(vip.HeatTemp_Vip_6[i])
    hvals = filter(!isnan, [h1, h2, h3, h5, h6])
    T_heat_C = isempty(hvals) ? 20.0 : mean(hvals)   # heating setpoint [°C]
    U11 = T_heat_C + 273.15

    # [D] Tpipe: measured if valid; else estimate from heating setpoint
    pg_c = (!isnan(pg) && pg > PIPE_MIN_VALID) ? pg : NaN
    pl_c = (!isnan(pl) && pl > PIPE_MIN_VALID) ? pl : NaN
    pipe_max = max(isnan(pg_c) ? -Inf : pg_c, isnan(pl_c) ? -Inf : pl_c)

    if pipe_max == -Inf
        To_val = isnan(To) ? 15.0 : To
        if To_val < T_heat_C - 2.0 && To_val < PIPE_WARM_TOUT
            # outdoor cold enough that heating circuit was likely running
            Tpipe = (T_heat_C + PIPE_DELTA_K) + 273.15
            global n_pipe_heated += 1
        else
            # warm outdoor: heating not needed, pipe at ambient
            Tpipe = To_val + 2.0 + 273.15
            global n_pipe_cold += 1
        end
    else
        Tpipe = pipe_max + 273.15
        global n_pipe_meas += 1
    end

    # [A]+[E] CO2 dosing: proportional fraction when setpoint above threshold
    co2vip = _num(vip.CO2_Vip[i])
    U10 = (isnan(co2vip) || co2vip <= CO2_AMB_PPM + CO2_THR_PPM) ? 0.0 : CO2_INJECT_FRAC

    # ── observations ────────────────────────────────────────────────────────
    Ta = _num(gh.Tair[i]);  Co = _num(gh.CO2air[i]);  Rh = _num(gh.RHair[i])

    push!(ts,       serial2unix(g))
    push!(Ig_v,     Igv);   push!(Tout_v, ToK);  push!(Wind_v, Wsv)
    push!(VPout_v,  _vpout(isnan(To) ? 15.0 : To, RHov))
    push!(Tsky_v,   Tskv)
    push!(U1v,  U1);    push!(U8v,  U8);    push!(U10v,  U10)
    push!(U11v, U11);   push!(U12v, U12);   push!(Tpipev, Tpipe)
    push!(Tair_obs, isnan(Ta) ? NaN : Ta + 273.15)
    push!(CO2_obs,  isnan(Co) ? NaN : Co * 536.4 / ((isnan(Ta) ? 20.0 : Ta) + 273.15))
    push!(RH_obs,   isnan(Rh) ? NaN : Rh)
end

@assert length(ts) > 100 "Too few rows found — check DAY0_SERIAL"

if median(Wind_v) > 25.0
    Wind_v ./= 3.6
    @info "Wind converted from km/h to m/s"
end

t0   = ts[1]
trel = ts .- t0

A(v) = stepfun(ts,   v)
Z(v) = stepfun(trel, v)
zc   = Z(zeros(length(trel)))

# ── Tpipe and dosing diagnostics ─────────────────────────────────────────────
ntot = n_pipe_meas + n_pipe_heated + n_pipe_cold
@printf("Tpipe source:  measured=%d (%.0f%%)  heat-est=%d (%.0f%%)  cold-fallback=%d (%.0f%%)\n",
    n_pipe_meas,    100n_pipe_meas/ntot,
    n_pipe_heated,  100n_pipe_heated/ntot,
    n_pipe_cold,    100n_pipe_cold/ntot)
frac_dosing = mean(U10v .> 0)
@printf("CO2 dosing ON: %.1f%% of timesteps  (U10 = %.3f when ON)\n",
    100frac_dosing, CO2_INJECT_FRAC)

# ── build SimBase ─────────────────────────────────────────────────────────────
params = load_params(CLM_JS, CRP_JS)

weather = WeatherInputs(
    A(Ig_v),   A(Tsky_v),  A(Tout_v),
    _t -> 273.15, _t -> 18.0 + 273.15,
    A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v),
    _t -> LAI_FIX)

controls = ControlInputs(
    Z(U1v),  zc,   zc,  zc,   zc,   zc,
    zc,      Z(U8v), zc, Z(U10v), Z(U11v), Z(U12v),
    Z(Tpipev))

base = SimBase(params, weather, controls, t0)

# ── initial state ─────────────────────────────────────────────────────────────
i0   = findfirst(!isnan, Tair_obs)
Ta0  = isnothing(i0) ? 296.15 : Tair_obs[i0]
Co0  = isnothing(i0) ? 575.0  : (isnan(CO2_obs[i0]) ? 575.0 : CO2_obs[i0])
Rh0  = isnothing(i0) ? 80.0   : (isnan(RH_obs[i0])  ? 80.0  : RH_obs[i0])
V1_0 = (Rh0 / 100.0) * Pws(Ta0)
u0   = [Ta0, Ta0, Co0, V1_0, 0.0, 18.0 + 273.15]

# ── simulate ──────────────────────────────────────────────────────────────────
days  = collect(0.0:(N_DAYS - 1))
tspan = (0.0, 86400.0)
teval = 0.0:300.0:86400.0

println("Running Vanthoor forward simulation for $N_DAYS days (v3: Tpipe + CO2 rate) …")
@time _, T2_sim, RH_sim, CO2_sim = forward_map(
    CAL_X, CAL_NAMES, base, days, tspan, teval; u0 = copy(u0))

# ── map observations to prediction grid ──────────────────────────────────────
predtimes = vcat(t0, [t0 + d*86400.0 + tj for d in days for tj in teval])

Aobs_T   = stepfun(ts, Tair_obs)
Aobs_CO2 = stepfun(ts, CO2_obs)
Aobs_RH  = stepfun(ts, RH_obs)

Tobs_g   = Aobs_T.(predtimes)
CO2obs_g = Aobs_CO2.(predtimes)
RHobs_g  = Aobs_RH.(predtimes)

T2v   = collect(T2_sim)
CO2v  = collect(CO2_sim)
RHv   = collect(RH_sim)

# ── overall RMSE ─────────────────────────────────────────────────────────────
mask_T  = .!isnan.(Tobs_g)
mask_C  = .!isnan.(CO2obs_g)
mask_RH = .!isnan.(RHobs_g)

rmse_T   = sqrt(mean((T2v[mask_T]   .- Tobs_g[mask_T]).^2))
rmse_CO2 = sqrt(mean((CO2v[mask_C]  .- CO2obs_g[mask_C]).^2))
rmse_RH  = sqrt(mean((RHv[mask_RH]  .- RHobs_g[mask_RH]).^2))
co2_ppm  = rmse_CO2 * (20.0+273.15) / 536.4

bias_T   = mean(T2v[mask_T] .- Tobs_g[mask_T])

println()
println("=" ^ 65)
println("  Vanthoor RMSE (v3: Tpipe-heat + CO2 rate) — days 0–$(N_DAYS-1) AiCU")
println("=" ^ 65)
@printf("  RMSE T_air  : %5.2f K       (v2: 4.80 K)\n",   rmse_T)
@printf("  Bias T_air  : %+6.2f K\n",   bias_T)
@printf("  RMSE CO2    : %5.1f mg/m³   (v2: 749.0)  (~%.0f ppm)\n", rmse_CO2, co2_ppm)
@printf("  RMSE RH     : %5.1f %%      (v2: 13.9 %%)\n",   rmse_RH)
println()

# ── [C] per-day breakdown ─────────────────────────────────────────────────────
n_per_day = length(teval) + 1
println(@sprintf("  %-6s  %-8s  %-8s  %-8s  %-8s  %-10s",
        "Day", "BiasT(K)", "RMSE_T", "MaxPred", "MeanU8", "Tpipe_est"))
println("  " * "-"^62)

for d in 0:(N_DAYS-1)
    idx = (d * n_per_day + 2):((d+1)*n_per_day + 1)
    idx = idx[idx .<= length(T2v)]
    isempty(idx) && continue
    T2d = T2v[idx]
    Tod = Tobs_g[idx]
    valid = .!isnan.(Tod)
    isempty(findall(valid)) && continue

    bias   = mean(T2d[valid] .- Tod[valid])
    rmse_d = sqrt(mean((T2d[valid] .- Tod[valid]).^2))
    maxT   = maximum(T2d) - 273.15

    day_start = t0 + d * 86400.0
    day_end   = day_start + 86400.0
    imask = (ts .>= day_start) .& (ts .< day_end)
    u8_d   = isempty(U8v[imask])    ? NaN : mean(U8v[imask])
    tpipe_d = isempty(Tpipev[imask]) ? NaN : mean(Tpipev[imask]) - 273.15

    @printf("  %-6d  %+7.2f   %7.2f   %7.1f°C  %6.3f   %6.1f°C\n",
            d, bias, rmse_d, maxT, u8_d, tpipe_d)
end

println()
if rmse_T < 1.0
    println("  ✓  RMSE_T < 1 K  →  Sequential inference adequate.")
elseif rmse_T < 2.0
    println("  ~  1–2 K  →  Borderline; calibrate nu4/gamma4 before joint inference.")
elseif rmse_T < 3.0
    println("  !  2–3 K  →  Run run_inference_venting.jl to calibrate parameters.")
elseif rmse_T < 4.0
    println("  !  3–4 K  →  Parameter calibration needed; check CO2_INJECT_FRAC tuning.")
else
    println("  ✗  RMSE_T > 4 K  →  Structural issue remains; revisit PIPE_DELTA_K or venting.")
end
println("=" ^ 65)
@printf("\nKey constants to tune if needed:\n")
@printf("  PIPE_DELTA_K    = %.1f K   (pipe above setpoint when heating; try 10–25)\n", PIPE_DELTA_K)
@printf("  CO2_INJECT_FRAC = %.3f    (CO2 injection rate fraction; try 0.01–0.20)\n", CO2_INJECT_FRAC)

# ── plots ─────────────────────────────────────────────────────────────────────
hours = (predtimes .- t0) ./ 3600.0

pT = plot(hours, Tobs_g .- 273.15;
          label="measured Tair", lc=:black, lw=2,
          xlabel="Hours since day 0", ylabel="T air (°C)",
          title="Vanthoor T_air v3 — RMSE=$(round(rmse_T,digits=2)) K, bias=$(round(bias_T,digits=2)) K",
          legend=:topleft, size=(900,350))
plot!(pT, hours, T2v .- 273.15; label="model v3", lc=:darkorange, lw=1.5)
savefig(pT, joinpath(FIGDIR, "rmse_v3_Tair.png"))

pC = plot(hours, CO2obs_g .* (20.0+273.15)./536.4;
          label="measured CO2", lc=:black, lw=2,
          xlabel="Hours since day 0", ylabel="CO₂ (ppm)",
          title="Vanthoor CO₂ v3 — RMSE=$(round(co2_ppm,digits=0)) ppm",
          legend=:topleft, size=(900,350))
plot!(pC, hours, CO2v .* (20.0+273.15)./536.4; label="model v3", lc=:steelblue, lw=1.5)
savefig(pC, joinpath(FIGDIR, "rmse_v3_CO2.png"))

pPipe = plot(hours[1:length(ts)], Tpipev .- 273.15;
             label="Tpipe used", lc=:firebrick, lw=1,
             xlabel="Hours since day 0", ylabel="T pipe (°C)",
             title="Pipe temperature input (v3)",
             legend=:topright, size=(900,250))
plot!(pPipe, hours[1:length(ts)], Tout_v .- 273.15; label="Tout", lc=:steelblue, lw=1, ls=:dash)
savefig(pPipe, joinpath(FIGDIR, "rmse_v3_Tpipe.png"))

println("\nFigures: rmse_v3_Tair.png, rmse_v3_CO2.png, rmse_v3_Tpipe.png")
