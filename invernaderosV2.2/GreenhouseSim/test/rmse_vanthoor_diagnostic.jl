# =============================================================================
# rmse_vanthoor_diagnostic.jl
#
# Run the V2.2 Vanthoor ODE at literature parameters on days 0–41
# (AiCU AGC-2018), driven by measured actuator states from Greenhouse_climate.csv
# and outdoor weather from meteo.csv.
# Compare T2 (air temperature) and CO2 to measured Tair and CO2air.
#
# RMSE < 1 K  → sequential inference (Note 25 Strategy 1) is adequate
# RMSE > 2 K  → joint inference (Strategy 2 Option B) required
#              (but first check vent-angle convention and run run_inference_venting.jl)
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
const DAY0_SERIAL = 43326.0     # Excel serial for 2018-08-14 (season day 0)
const N_DAYS      = 42          # days 0–41
const LAI_FIX     = 2.0         # fixed LAI for climate run
const EXCEL_UNIX0 = 25569.0     # Excel epoch offset
const SIGMA_SB    = 5.670374419e-8

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
U1v      = Float64[];  U8v    = Float64[];  U11v   = Float64[]
U12v     = Float64[];  Tpipev = Float64[]
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
    # sky temperature from net longwave (Pyrgeo)
    Ldn  = SIGMA_SB * ToK^4 + (isnan(Py) ? 0.0 : Py)
    Tskv = isnan(Py) ? (-0.4 + 273.15) : clamp((Ldn / SIGMA_SB)^0.25, ToK-30.0, ToK+2.0)

    # ── measured actuator states from GH_climate ─────────────────────────────
    vl = _num(gh.VentLee[i]);   vw  = _num(gh.Ventwind[i])
    en = _num(gh.EnScr[i]);     al  = _num(gh.AssimLight[i])
    pg = _num(gh.PipeGrow[i]);  pl  = _num(gh.PipeLow[i])

    U8  = clamp(((isnan(vl) ? 0.0 : vl) + (isnan(vw) ? 0.0 : vw)) / 2.0 / 100.0, 0.0, 1.0)
    U1  = isnan(en) ? 0.0 : clamp(en  / 100.0, 0.0, 1.0)
    U12 = isnan(al) ? 0.0 : clamp(al  / 100.0, 0.0, 1.0)
    pipe_max = max(isnan(pg) ? -Inf : pg, isnan(pl) ? -Inf : pl)
    Tpipe = (pipe_max == -Inf) ? ToK : (pipe_max + 273.15)

    # ── U11: heating setpoint from vip.csv (mean of active zones 1,2,3,5,6) ─
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
    push!(U1v,  U1);    push!(U8v,  U8);    push!(U11v, U11)
    push!(U12v, U12);   push!(Tpipev, Tpipe)
    push!(Tair_obs, isnan(Ta) ? NaN : Ta + 273.15)
    push!(CO2_obs,  isnan(Co) ? NaN : Co * 536.4 / ((isnan(Ta) ? 20.0 : Ta) + 273.15))
    push!(RH_obs,   isnan(Rh) ? NaN : Rh)
end

@assert length(ts) > 100 "Too few rows found in window — check DAY0_SERIAL"

# fix wind units: meteo.csv may be km/h
if median(Wind_v) > 25.0
    Wind_v ./= 3.6
    @info "Wind converted from km/h to m/s"
end

t0   = ts[1]
trel = ts .- t0

A(v) = stepfun(ts,   v)    # indexed by absolute Unix time
Z(v) = stepfun(trel, v)    # indexed by seconds since experiment start
zc   = Z(zeros(length(trel)))

# ── build SimBase ─────────────────────────────────────────────────────────────
params = load_params(CLM_JS, CRP_JS)

weather = WeatherInputs(
    A(Ig_v),   A(Tsky_v),  A(Tout_v),
    _t -> 273.15, _t -> 18.0 + 273.15,
    A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v),
    _t -> LAI_FIX)

controls = ControlInputs(
    Z(U1v),  zc,   zc,  zc,   zc,   zc,
    zc,      Z(U8v), zc, zc, Z(U11v), Z(U12v),
    Z(Tpipev))

base = SimBase(params, weather, controls, t0)

# ── initial state from first valid observation ────────────────────────────────
i0   = findfirst(!isnan, Tair_obs)
Ta0  = isnothing(i0) ? 296.15 : Tair_obs[i0]
Co0  = isnothing(i0) ? 575.0  : (isnan(CO2_obs[i0]) ? 575.0 : CO2_obs[i0])
Rh0  = isnothing(i0) ? 80.0   : (isnan(RH_obs[i0])  ? 80.0  : RH_obs[i0])
V1_0 = (Rh0 / 100.0) * Pws(Ta0)
u0   = [Ta0, Ta0, Co0, V1_0, 0.0, 18.0 + 273.15]

# ── simulate ──────────────────────────────────────────────────────────────────
days  = collect(0.0:(N_DAYS - 1))
tspan = (0.0, 86400.0)
teval = 0.0:300.0:86400.0     # 5-min output, matches observations

println("Running Vanthoor forward simulation for $N_DAYS days …")
@time _, T2_sim, RH_sim, CO2_sim = forward_map(
    Float64[], String[], base, days, tspan, teval; u0 = copy(u0))

# ── interpolate observations onto prediction time grid ───────────────────────
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

# ── RMSE ─────────────────────────────────────────────────────────────────────
mask_T   = .!isnan.(Tobs_g)
mask_C   = .!isnan.(CO2obs_g)
mask_RH  = .!isnan.(RHobs_g)

rmse_T   = sqrt(mean((T2v[mask_T]   .- Tobs_g[mask_T]).^2))
rmse_CO2 = sqrt(mean((CO2v[mask_C]  .- CO2obs_g[mask_C]).^2))
rmse_RH  = sqrt(mean((RHv[mask_RH]  .- RHobs_g[mask_RH]).^2))
co2_ppm  = rmse_CO2 * (20.0 + 273.15) / 536.4   # rough ppm conversion at 20°C

println()
println("=" ^ 60)
println("  Vanthoor RMSE diagnostic  —  days 0–$(N_DAYS-1) (AiCU, literature params)")
println("=" ^ 60)
@printf("  RMSE T_air  : %5.2f K       (threshold: ~1 K for sequential)\n", rmse_T)
@printf("  RMSE CO2    : %5.1f mg/m³   (~%.0f ppm at 20°C)\n", rmse_CO2, co2_ppm)
@printf("  RMSE RH     : %5.1f %%\n",   rmse_RH)
println()
if rmse_T < 1.0
    println("  ✓  RMSE_T < 1 K  →  Sequential inference (Note 25 Strategy 1) is adequate.")
elseif rmse_T < 2.0
    println("  ~  1 K < RMSE_T < 2 K  →  Borderline. Check vent-angle convention,")
    println("     then run run_inference_venting.jl to update nu4/gamma4.")
else
    println("  ✗  RMSE_T > 2 K  →  Joint inference likely needed (Note 25 Strategy 2 Option B).")
    println("     First check: is VentLee divided by 100 correct? (full-open angle in degrees)")
end
println("=" ^ 60)

# ── plots ─────────────────────────────────────────────────────────────────────
hours = (predtimes .- t0) ./ 3600.0

# Air temperature
pT = plot(hours, Tobs_g .- 273.15;
          label = "measured Tair", lc = :black, lw = 2,
          xlabel = "Hours since day 0", ylabel = "T air (°C)",
          title  = "Vanthoor T_air — days 0–$(N_DAYS-1)  (RMSE = $(round(rmse_T,digits=2)) K)",
          legend = :topleft, size = (900, 350))
plot!(pT, hours, T2v .- 273.15; label = "model (literature params)", lc = :darkorange, lw = 1.5)
savefig(pT, joinpath(FIGDIR, "rmse_vanthoor_Tair.png"))

# CO2
pC = plot(hours, CO2obs_g .* (20.0+273.15)./536.4;
          label = "measured CO2", lc = :black, lw = 2,
          xlabel = "Hours since day 0", ylabel = "CO₂ (ppm)",
          title  = "Vanthoor CO₂ — RMSE = $(round(co2_ppm,digits=0)) ppm",
          legend = :topleft, size = (900, 350))
plot!(pC, hours, CO2v .* (20.0+273.15)./536.4; label = "model", lc = :steelblue, lw = 1.5)
savefig(pC, joinpath(FIGDIR, "rmse_vanthoor_CO2.png"))

println("\nFigures saved to: $(FIGDIR)")
println("  rmse_vanthoor_Tair.png")
println("  rmse_vanthoor_CO2.png")
