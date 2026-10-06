# =============================================================================
# rmse_vanthoor_diagnostic.jl  (v2 — CO2 dosing + pipe-floor clamping)
#
# Changes vs. v1:
#   [A] CO2 dosing: U10 driven by CO2_Vip setpoint (zone-averaged) when above
#       ambient (664 mg/m³ ≈ 390 ppm at 20°C).  psi2=13300 mg/m²/s is the
#       dosing valve capacity; U10 ∈ [0,1] is a fraction of that capacity.
#       Simple bang-bang: U10 = 1 if CO2_Vip > 440 ppm (≈390+50), else 0. CO2_Vip is in PPM.
#       (A PID would be more realistic but needs tuning; bang-bang explores
#       whether the CO2 channel matters at all for the temperature RMSE.)
#
#   [B] Pipe temperature clamp: PipeGrow/PipeLow = 0 means "heating off, pipe
#       at ambient", NOT 0 °C.  Clamp Tpipe = max(reading, Tout+2) K so pipes
#       never artificially cool the air.
#
#   [C] Daily-breakdown table: print per-day mean T_air bias and max T2.
#
# Run from the GreenhouseSim folder:
#   julia --project=. test/rmse_vanthoor_diagnostic.jl
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Statistics, Printf, Plots

# ── paths ────────────────────────────────────────────────────────────────────
const INV    = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TEAM   = "AiCU"
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

# CO2 dosing: dose when setpoint exceeds ambient by > threshold
const CO2_AMB_PPM = 390.0    # ≈outside CO₂ in ppm
const CO2_THR_PPM = 50.0     # start dosing when setpoint > ambient + threshold [ppm]

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
U1v      = Float64[];  U8v    = Float64[];  U10v   = Float64[]
U11v     = Float64[];  U12v   = Float64[];  Tpipev = Float64[]
Ig_v     = Float64[];  Tout_v = Float64[];  Wind_v = Float64[]
VPout_v  = Float64[];  Tsky_v = Float64[]
Tair_obs = Float64[];  CO2_obs = Float64[]; RH_obs = Float64[]

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

    # ── measured actuator states ─────────────────────────────────────────────
    vl = _num(gh.VentLee[i]);   vw  = _num(gh.Ventwind[i])
    en = _num(gh.EnScr[i]);     al  = _num(gh.AssimLight[i])
    pg = _num(gh.PipeGrow[i]);  pl  = _num(gh.PipeLow[i])

    U8  = clamp(((isnan(vl) ? 0.0 : vl) + (isnan(vw) ? 0.0 : vw)) / 2.0 / 100.0, 0.0, 1.0)
    U1  = isnan(en) ? 0.0 : clamp(en  / 100.0, 0.0, 1.0)
    U12 = isnan(al) ? 0.0 : clamp(al  / 100.0, 0.0, 1.0)

    # [B] Pipe clamp: readings ≤ 5°C mean "off" → use Tout+2 K (ambient-ish)
    pg_c = (!isnan(pg) && pg > 5.0) ? pg : NaN
    pl_c = (!isnan(pl) && pl > 5.0) ? pl : NaN
    pipe_max = max(isnan(pg_c) ? -Inf : pg_c, isnan(pl_c) ? -Inf : pl_c)
    Tpipe = (pipe_max == -Inf) ? (isnan(To) ? 15.0 : To) + 2.0 + 273.15 : pipe_max + 273.15

    # [A] CO2 dosing: U10=1 when CO2_Vip setpoint > ambient + threshold
    co2vip = _num(vip.CO2_Vip[i])
    U10 = (isnan(co2vip) || co2vip <= CO2_AMB_PPM + CO2_THR_PPM) ? 0.0 : 1.0

    # U11: heating setpoint from vip.csv (mean of zones 1,2,3,5,6)
    h1 = _num(vip.HeatTemp_Vip_1[i]); h2 = _num(vip.HeatTemp_Vip_2[i])
    h3 = _num(vip.HeatTemp_Vip_3[i]); h5 = _num(vip.HeatTemp_Vip_5[i])
    h6 = _num(vip.HeatTemp_Vip_6[i])
    hvals = filter(!isnan, [h1, h2, h3, h5, h6])
    U11 = (isempty(hvals) ? 20.0 : mean(hvals)) + 273.15

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

# ── CO2 dosing stats ─────────────────────────────────────────────────────────
frac_dosing = mean(U10v .> 0)
println("CO2 dosing ON: $(round(100*frac_dosing, digits=1))% of timesteps")

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

println("Running Vanthoor forward simulation for $N_DAYS days (v2: CO2 + pipe fix) …")
@time _, T2_sim, RH_sim, CO2_sim = forward_map(
    Float64[], String[], base, days, tspan, teval; u0 = copy(u0))

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

println()
println("=" ^ 65)
println("  Vanthoor RMSE (v2: CO2+pipe fix) — days 0–$(N_DAYS-1) AiCU")
println("=" ^ 65)
@printf("  RMSE T_air  : %5.2f K       (v1: 5.64 K)\n",   rmse_T)
@printf("  RMSE CO2    : %5.1f mg/m³   (v1: 729.1)  (~%.0f ppm)\n", rmse_CO2, co2_ppm)
@printf("  RMSE RH     : %5.1f %%      (v1: 13.2 %%)\n",   rmse_RH)
println()

# ── [C] per-day breakdown ─────────────────────────────────────────────────────
n_per_day = length(teval) + 1    # 5-min for 1 day
println(@sprintf("  %-6s  %-8s  %-8s  %-8s  %-8s  %-8s",
        "Day", "BiasT(K)", "RMSE_T", "MaxPred", "MeanU8", "MeanVent"))
println("  " * "-"^58)

for d in 0:(N_DAYS-1)
    # slice prediction
    idx = (d * n_per_day + 2):((d+1)*n_per_day + 1)
    idx = idx[idx .<= length(T2v)]
    isempty(idx) && continue
    T2d = T2v[idx]
    Tod = Tobs_g[idx]
    valid = .!isnan.(Tod)
    isempty(findall(valid)) && continue

    bias = mean(T2d[valid] .- Tod[valid])
    rmse_d = sqrt(mean((T2d[valid] .- Tod[valid]).^2))
    maxT   = maximum(T2d) - 273.15

    # mean U8 for this day from the forcing
    day_start = t0 + d * 86400.0
    day_end   = day_start + 86400.0
    imask = (ts .>= day_start) .& (ts .< day_end)
    u8_d = isempty(U8v[imask]) ? NaN : mean(U8v[imask])
    vl_d = isempty(U8v[imask]) ? NaN : mean(U8v[imask]) * 200   # approx VentLee%

    @printf("  %-6d  %+7.2f   %7.2f   %7.1f°C  %6.3f   %5.1f%%\n",
            d, bias, rmse_d, maxT, u8_d, vl_d)
end

println()
if rmse_T < 1.0
    println("  ✓  RMSE_T < 1 K  →  Sequential inference adequate (Note 25 Strategy 1).")
elseif rmse_T < 2.0
    println("  ~  1–2 K  →  Borderline; run run_inference_venting.jl to calibrate.")
elseif rmse_T < 4.0
    println("  !  2–4 K  →  Parameter calibration likely needed before joint inference.")
else
    println("  ✗  RMSE_T > 4 K  →  Structural issue; check venting / radiation params.")
end
println("=" ^ 65)

# ── plots ─────────────────────────────────────────────────────────────────────
hours = (predtimes .- t0) ./ 3600.0

pT = plot(hours, Tobs_g .- 273.15;
          label="measured Tair", lc=:black, lw=2,
          xlabel="Hours since day 0", ylabel="T air (°C)",
          title="Vanthoor T_air v2 — RMSE=$(round(rmse_T,digits=2)) K",
          legend=:topleft, size=(900,350))
plot!(pT, hours, T2v .- 273.15; label="model v2", lc=:darkorange, lw=1.5)
savefig(pT, joinpath(FIGDIR, "rmse_v2_Tair.png"))

pC = plot(hours, CO2obs_g .* (20.0+273.15)./536.4;
          label="measured CO2", lc=:black, lw=2,
          xlabel="Hours since day 0", ylabel="CO₂ (ppm)",
          title="Vanthoor CO₂ v2 — RMSE=$(round(co2_ppm,digits=0)) ppm",
          legend=:topleft, size=(900,350))
plot!(pC, CO2v .* (20.0+273.15)./536.4; label="model v2", lc=:steelblue, lw=1.5)
savefig(pC, joinpath(FIGDIR, "rmse_v2_CO2.png"))

println("\nFigures: rmse_v2_Tair.png, rmse_v2_CO2.png")
