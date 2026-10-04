# =============================================================================
# validate_leakage_agc.jl  --  does the CORRECTED model reproduce measured Tair
#                              on a real VENTING day WITHOUT inflated leakage?
# =============================================================================
# Run:  julia --project=. test/validate_leakage_agc.jl
#
# Background. In V2.2 climate.jl we fixed two ventilation bugs (the vent-rate sign
# and the missing roof-vent term f4 in the air heat balance h7). The heat leakage
# nu4 had been inflated from Vanthoor's physical ~1e-4 to 0.01 *purely to stabilize
# the ODE* (see the JSON comment) -- i.e. to paper over those bugs. The twin sweep
# already showed the twin is stable and better-behaved at 1e-4. This script is the
# DATA-based check: drive the corrected climate model with the Autonomous Greenhouse
# Challenge 2018 measured weather + measured ROOF-VENT apertures over a warm venting
# window, and compare predicted air temperature against the MEASURED air temperature
# at nu4 = 1e-4 (physical) vs 1e-3 vs 1e-2 (old inflated). If the physical value
# tracks the data as well as the inflated one, the inflation is unnecessary.
#
# Modelling choices (documented; adjust at the top if you know better values):
#   * VENTS: the AGC compartment is a Venlo house with two ROOF vents (leeward /
#     windward sides of the ridge) and no side walls. Both map to the model's ROOF
#     term U8 = mean(VentLee,Ventwind)/100; the side/wind vent U6 = 0. (The roof
#     term f7 already carries the wind contribution via its I8^2 term, so using U6
#     for a roof vent would double-count wind.) THIS is the term the fix corrects.
#   * PIPES: I3 = max(mean(PipeGrow,PipeLow), Tout) in K, so an "off" (~0) pipe acts
#     as no heat source rather than a cold sink. Minor on a venting day (pipes ~0).
#   * LAI: fixed at LAI_FIX for the window (transpiration cooling). Its exact value
#     shifts all curves together, so it does NOT bias the 1e-4-vs-1e-2 comparison.
#   * Tsky/Tsoil/CO2out: constants (daytime venting fit; FIR/soil/CO2 are second-order
#     and CO2 does not enter the T/RH balances at all).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
const GS = GreenhouseControl.GreenhouseSim      # reach engine internals directly
using Plots, Printf, Statistics

# ------------------------------------------------------------------ user knobs
INV      = normpath(joinpath(@__DIR__, "..", "..", ".."))
team_dir = joinpath(INV, "invernaderos2.0", "data", "Reference(Growers)")
gh_csv   = joinpath(team_dir, "Greenhouse_climate.csv")
meteo    = joinpath(INV, "invernaderos2.0", "data", "meteo.csv")
clim_js  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_js  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")

WIN_START_SERIAL = 43346.0     # Excel-serial day the window starts (day 20; venting days 21-22)
NWIN_DAYS        = 3           # window length (days)
LAI_FIX          = 2.3         # standing LAI ~season day 20-23 (calibrated crop)
NU4_LIST         = [1.0e-4, 1.0e-3, 1.0e-2]   # physical -> old inflated

EXCEL_UNIX0 = 25569.0          # Excel serial of 1970-01-01
serial2unix(s) = (s - EXCEL_UNIX0) * 86400.0
SIGMA_SB = 5.670374419e-8      # Stefan-Boltzmann, for sky temp from net longwave

# --------------------------------------------------------------- small helpers
_num(x) = (x isa Number && !(x isa Missing) && x == x) ? Float64(x) : NaN

"Step-function closure over sorted knot times `ts` (constant extrapolation)."
function stepfun(ts::Vector{Float64}, vs::Vector{Float64})
    return t -> @inbounds vs[clamp(searchsortedlast(ts, t), 1, length(ts))]
end

using CSV, DataFrames

# --------------------------------------------------- read window rows (aligned)
gh = DataFrame(CSV.File(gh_csv))
mt = DataFrame(CSV.File(meteo))
n  = min(nrow(gh), nrow(mt))

win_lo = WIN_START_SERIAL
win_hi = WIN_START_SERIAL + NWIN_DAYS

ts=Float64[]; Tout=Float64[]; Ig=Float64[]; Wind=Float64[]; VPout=Float64[]; Tsky=Float64[]
U8=Float64[]; U1=Float64[]; U12=Float64[]; Tpipe=Float64[]
Tair_obs=Float64[]; RH_obs=Float64[]; CO2_obs=Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]);  (isnan(g) || g < win_lo || g >= win_hi) && continue
    Ta = _num(gh.Tair[i]);   isnan(Ta) && continue
    Rh = _num(gh.RHair[i]);  Co = _num(gh.CO2air[i])
    vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i])
    pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    en = _num(gh.EnScr[i]);   al = _num(gh.AssimLight[i])
    To = _num(mt.Tout[i]);    Igv = _num(mt.Iglob[i]); Ws = _num(mt.Windsp[i]); Rho = _num(mt.Rhout[i])
    Py = _num(mt.Pyrgeo[i])   # net longwave [W/m2] (negative = radiative loss)

    push!(ts, serial2unix(g))
    ToK = (isnan(To) ? 15.0 : To) + 273.15
    push!(Tout, ToK)
    push!(Ig,   isnan(Igv) ? 0.0 : max(Igv, 0.0))
    push!(Wind, isnan(Ws) ? 2.0 : Ws)
    ToC = isnan(To) ? 15.0 : To
    rh  = isnan(Rho) ? 80.0 : Rho
    push!(VPout, rh/100 * GS.Pws(ToC + 273.15))
    # sky temperature from MEASURED net longwave: L_down = sigma*Tout^4 + Pyrgeo,
    # Tsky = (L_down/sigma)^0.25. Warm (~Tout) on cloudy nights, cold on clear ones.
    Ldn  = SIGMA_SB * ToK^4 + (isnan(Py) ? 0.0 : Py)
    Tsk  = Ldn > 0 ? (Ldn / SIGMA_SB)^0.25 : (ToK - 20.0)
    push!(Tsky, isnan(Py) ? (-0.4 + 273.15) : clamp(Tsk, ToK - 30.0, ToK + 2.0))
    vlv = isnan(vl) ? 0.0 : vl; vwv = isnan(vw) ? 0.0 : vw
    push!(U8, clamp((vlv + vwv)/2 / 100, 0.0, 1.0))       # <-- roof aperture fraction
    push!(U1, isnan(en) ? 0.0 : clamp(en/100, 0.0, 1.0))  # thermal screen
    push!(U12, isnan(al) ? 0.0 : clamp(al/100, 0.0, 1.0)) # lamps
    pipe = maximum(skipmissing([isnan(pg) ? -Inf : pg, isnan(pl) ? -Inf : pl]))
    pipeK = (pipe == -Inf) ? (isnan(To) ? 288.15 : To+273.15) : pipe + 273.15
    push!(Tpipe, max(pipeK, (isnan(To) ? 288.15 : To+273.15)))  # off pipe -> no cold sink
    push!(Tair_obs, Ta + 273.15)
    push!(RH_obs, isnan(Rh) ? NaN : Rh)
    push!(CO2_obs, isnan(Co) ? 575.0 : Co * 536.4 / (Ta + 273.15))
end
@assert length(ts) > 10 "no rows found in window; check WIN_START_SERIAL"

# wind units guard: AGC Windsp should be m/s; if it looks like km/h, convert
if median(Wind) > 25; Wind ./= 3.6; @info "windspeed looked like km/h; converted to m/s"; end

t0 = ts[1]
trel = ts .- t0                       # controls: seconds since window start
Z(v) = stepfun(trel, v)               # control interpolants (relative time)
A(v) = stepfun(ts, v)                 # weather interpolants (absolute unix)

# ------------------------------------------------------------- build SimBase
params = load_params(clim_js, crop_js)

weather = GS.WeatherInputs(
    A(Ig),                    # Iglobal
    A(Tsky),                  # Tsky (from measured net longwave, per-hour)
    A(Tout),                  # Tout
    _t -> 273.15,             # TmechCool (unused)
    _t -> 18.0 + 273.15,      # Tsoil
    A(Wind),                  # WindSpeed
    _t -> 834.7,              # CO2out (does not affect T/RH)
    A(VPout),                 # VPout
    A(Ig),                    # Idocel
    _t -> LAI_FIX,            # LAI
)

zero_c = Z(zeros(length(trel)))
controls = GS.ControlInputs(
    Z(U1),  zero_c, zero_c, zero_c, zero_c, zero_c,   # U1(screen) U2..U5 U6(side)=0
    zero_c, Z(U8),  zero_c, zero_c, zero_c, Z(U12),   # U7 U8(roof) U9 U10 U11 U12(lamps)
    Z(Tpipe),
)

# initial state from first measured sample
V1_0 = (isnan(RH_obs[1]) ? 80.0 : RH_obs[1])/100 * GS.Pws(Tair_obs[1])
u0   = [Tair_obs[1], Tair_obs[1], CO2_obs[1], V1_0, 0.0, 18.0+273.15]  # last = floor T5

days  = collect(0.0:(NWIN_DAYS-1))
tspan = (0.0, 86400.0)
teval = 0.0:300.0:86400.0

# reconstruct absolute timestamps of the predicted grid (u0, then each day's teval)
predtimes = vcat(t0, [t0 + d*86400.0 + tj for d in days for tj in teval])
Tobs_i = A(Tair_obs); RHobs_i = A(RH_obs)
Tmeas  = Tobs_i.(predtimes) .- 273.15
RHmeas = RHobs_i.(predtimes)
hours  = (predtimes .- t0) ./ 3600.0

# -------------------------------------------------------------- run + compare
base = GS.SimBase(params, weather, controls, t0)
pT = plot(title = "Air temperature: model vs measured (venting window)",
          xlabel = "hours from window start", ylabel = "T air (C)", legend = :topright)
plot!(pT, hours, Tmeas, lw = 3, lc = :black, label = "measured Tair")

println(@sprintf("%-9s %10s %12s", "nu4", "RMSE_all", "RMSE_daytime"))
daymask = [11.0 <= mod(h,24) <= 18.0 for h in hours]   # daytime, venting hours
for lk in NU4_LIST
    T1,T2,RH,CO2 = GS.forward_map([lk,lk], ["nu4","nu4CO2"], base, days, tspan, teval; u0 = copy(u0))
    Tp = collect(T2) .- 273.15
    m  = .!isnan.(Tmeas)
    rmse_all = sqrt(mean((Tp[m] .- Tmeas[m]).^2))
    dm = m .& daymask
    rmse_day = sqrt(mean((Tp[dm] .- Tmeas[dm]).^2))
    println(@sprintf("%-9.1e %10.2f %12.2f", lk, rmse_all, rmse_day))
    plot!(pT, hours, Tp, lw = 2, label = @sprintf("model nu4=%.0e", lk))
end

figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
savefig(pT, joinpath(figdir, "validate_leakage_agc.png"))
println("\nFigure: ", joinpath(figdir, "validate_leakage_agc.png"))
println("""
Reading it: if the nu4=1e-4 curve tracks the measured black line about as well as
(or better than) nu4=1e-2 -- especially through the daytime venting hours -- then the
inflated leakage is unnecessary under the corrected physics, and 1e-4 (Vanthoor) is
justified by the data. That closes the loop opened by the roof-vent fix.
""")
