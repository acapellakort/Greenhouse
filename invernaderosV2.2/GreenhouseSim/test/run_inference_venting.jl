# =============================================================================
# run_inference_venting.jl -- t-walk MCMC for nu4 (heat leakage) + gamma4 on a
#                             REAL venting day, under the corrected V2.2 physics
# =============================================================================
# Run from the GreenhouseSim folder:
#     julia --project=. test/run_inference_venting.jl
#
# Why this exists. The original run_inference.jl fits a cool, roof-closed window
# (U8 = 0 throughout), so it is inert to the roof-vent fix and cannot constrain
# leakage against ventilation. This script fits a day when the roof vents are wide
# open (Autonomous Greenhouse Challenge 2018, Reference/Growers), so the corrected
# air-heat balance h7 (which now includes the roof term f4) is actually exercised,
# and the heat-leakage nu4 is estimated where it matters.
#
# Design choices:
#   * VENTS: Venlo house, two ROOF vents -> U8 = mean(VentLee,Ventwind)/100, side
#     vent U6 = 0 (the roof term f7 already carries wind via its I8^2 term).
#   * Tsky from the MEASURED net longwave (Pyrgeo): realistic radiative loss.
#   * Likelihood scored ONLY over daytime venting hours (DAY_LO..DAY_HI). The ODE
#     still integrates the whole day, but night points are excluded so the model's
#     (thermal-mass) night bias cannot leak into the posterior.
#   * Observables: air temperature T2 and RH (RH pins gamma4). No canopy T (not
#     measured) and no CO2 (no dosing actuator logged; CO2 does not enter T/RH).
#
# Edit WIN_START_SERIAL / NWIN_DAYS / QoI_dict to change the day or parameter set.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, JTwalk, Plots, Printf, Statistics

# ------------------------------------------------------------------ user knobs
INV      = normpath(joinpath(@__DIR__, "..", "..", ".."))
team_dir = joinpath(INV, "invernaderos2.0", "data", "Reference(Growers)")
gh_csv   = joinpath(team_dir, "Greenhouse_climate.csv")
meteo    = joinpath(INV, "invernaderos2.0", "data", "meteo.csv")
clim_js  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_js  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")

WIN_START_SERIAL = 43347.0      # best single venting day (VentLee ~85-100 all day)
NWIN_DAYS        = 1
LAI_FIX          = 2.3
DAY_LO, DAY_HI   = 7.0, 20.0    # daytime venting hours scored in the likelihood
SIG_T, SIG_RH    = 1.5, 8.0     # observation std [C], [%RH]

chain_length = 120000
burn_in      = 40000

QoI = ["nu4", "gamma4", "k_ground"]
QoI_dict = Dict(                          # [min, prior-center, max]
    "nu4"    => [1.0e-5, 1.0e-4, 1.0e-1],  # heat leakage (Vanthoor ~1e-4; old inflated 1e-2)
    "gamma4" => [20.0,   82.0,   300.0],   # min stomatal resistance
    "k_ground"    => [0.5,   5.0,   50.0],     # single air<->ground conductance (fixed capacity C_soil)
)

EXCEL_UNIX0 = 25569.0
serial2unix(s) = (s - EXCEL_UNIX0) * 86400.0
SIGMA_SB = 5.670374419e-8
_num(x) = (x isa Number && !(x isa Missing) && x == x) ? Float64(x) : NaN
stepfun(ts, vs) = t -> @inbounds vs[clamp(searchsortedlast(ts, t), 1, length(ts))]

# --------------------------------------------------- read the venting window
gh = DataFrame(CSV.File(gh_csv)); mt = DataFrame(CSV.File(meteo))
n  = min(nrow(gh), nrow(mt))
win_lo = WIN_START_SERIAL; win_hi = WIN_START_SERIAL + NWIN_DAYS

ts=Float64[]; Tout=Float64[]; Ig=Float64[]; Wind=Float64[]; VPout=Float64[]; Tsky=Float64[]
U8=Float64[]; U1=Float64[]; U12=Float64[]; Tpipe=Float64[]
Tair_obs=Float64[]; RH_obs=Float64[]; CO2_obs=Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    Ta = _num(gh.Tair[i]);  isnan(Ta) && continue
    Rh = _num(gh.RHair[i]);  Co = _num(gh.CO2air[i])
    vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i])
    pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    en = _num(gh.EnScr[i]);   al = _num(gh.AssimLight[i])
    To = _num(mt.Tout[i]); Igv = _num(mt.Iglob[i]); Ws = _num(mt.Windsp[i]); Rho = _num(mt.Rhout[i]); Py = _num(mt.Pyrgeo[i])

    push!(ts, serial2unix(g))
    ToK = (isnan(To) ? 15.0 : To) + 273.15; push!(Tout, ToK)
    push!(Ig, isnan(Igv) ? 0.0 : max(Igv, 0.0))
    push!(Wind, isnan(Ws) ? 2.0 : Ws)
    push!(VPout, (isnan(Rho) ? 80.0 : Rho)/100 * Pws(ToK))
    Ldn = SIGMA_SB * ToK^4 + (isnan(Py) ? 0.0 : Py)
    Tsk = Ldn > 0 ? (Ldn/SIGMA_SB)^0.25 : ToK - 20.0
    push!(Tsky, isnan(Py) ? (-0.4+273.15) : clamp(Tsk, ToK-30.0, ToK+2.0))
    vlv = isnan(vl) ? 0.0 : vl; vwv = isnan(vw) ? 0.0 : vw
    push!(U8, clamp((vlv+vwv)/2/100, 0.0, 1.0))
    push!(U1, isnan(en) ? 0.0 : clamp(en/100, 0.0, 1.0))
    push!(U12, isnan(al) ? 0.0 : clamp(al/100, 0.0, 1.0))
    pipe = max(isnan(pg) ? -Inf : pg, isnan(pl) ? -Inf : pl)
    pipeK = pipe == -Inf ? ToK : pipe + 273.15
    push!(Tpipe, max(pipeK, ToK))
    push!(Tair_obs, Ta + 273.15)
    push!(RH_obs, isnan(Rh) ? NaN : Rh)
    push!(CO2_obs, isnan(Co) ? 575.0 : Co*536.4/(Ta+273.15))
end
@assert length(ts) > 20 "no rows in window"
if median(Wind) > 25; Wind ./= 3.6; end

t0 = ts[1]; trel = ts .- t0
Z(v) = stepfun(trel, v); A(v) = stepfun(ts, v)

params = load_params(clim_js, crop_js)
weather = WeatherInputs(A(Ig), A(Tsky), A(Tout), _t->273.15, _t->18.0+273.15,
                        A(Wind), _t->834.7, A(VPout), A(Ig), _t->LAI_FIX)
zc = Z(zeros(length(trel)))
controls = ControlInputs(Z(U1), zc, zc, zc, zc, zc, zc, Z(U8), zc, zc, zc, Z(U12), Z(Tpipe))
base = SimBase(params, weather, controls, t0)

# prediction grid + measured series sampled on it; daytime venting mask
days = collect(0.0:(NWIN_DAYS-1)); tspan = (0.0, 86400.0); teval = 0.0:900.0:86400.0
predtimes = vcat(t0, [t0 + d*86400.0 + tj for d in days for tj in teval])
hours = mod.((predtimes .- t0)./3600.0, 24.0)
Tmeas = A(Tair_obs).(predtimes); RHmeas = A(RH_obs).(predtimes)
mask = (hours .>= DAY_LO) .& (hours .<= DAY_HI) .& .!isnan.(Tmeas) .& .!isnan.(RHmeas)
@printf("scored points (daytime venting): %d of %d\n", sum(mask), length(mask))

V1_0 = (isnan(RH_obs[1]) ? 80.0 : RH_obs[1])/100 * Pws(Tair_obs[1])
u0   = [Tair_obs[1], Tair_obs[1], CO2_obs[1], V1_0, 0.0, 18.0+273.15]  # last = floor T5

# ------------------------------------------------------------------- MCMC
PriorSupp(x) = all(QoI_dict[QoI[k]][1] < x[k] < QoI_dict[QoI[k]][3] for k in eachindex(x))
function energy(x)
    _, T2, RH, _ = forward_map(x, QoI, base, days, tspan, teval; u0 = copy(u0))
    Tp = collect(T2); RHp = collect(RH)
    logL = -sum(((Tmeas[mask] .- Tp[mask])./SIG_T).^2) -
            sum(((RHmeas[mask] .- RHp[mask])./SIG_RH).^2)
    return -logL
end

n_p = length(QoI)
obj = jtwalk(n = n_p, U = energy, Supp = PriorSupp)
x0  = [rand()*(QoI_dict[q][3]-QoI_dict[q][1])+QoI_dict[q][1] for q in QoI]
xp0 = [rand()*(QoI_dict[q][3]-QoI_dict[q][1])+QoI_dict[q][1] for q in QoI]
println("Inferring: ", QoI, "  chain=", chain_length)
@time Run!(obj, T = chain_length, x0 = x0, xp0 = xp0)

# ------------------------------------------------------------------- report
figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
etr = obj.Output[:, end]
pE = plot(etr, title="Energy U=-logL", xlabel="iter", ylabel="U", label="")
vline!(pE, [burn_in], lc=:red, ls=:dash, label="burn-in")
savefig(pE, joinpath(figdir, "venting_energy.png"))

println("\n--- posteriors (post burn-in) ---")
med = Float64[]
for (i,name) in enumerate(QoI)
    c = obj.Output[burn_in:end, i]
    q = quantile(c, [0.025, 0.5, 0.975]); push!(med, q[2])
    @printf("%-8s  median=%.3e   95%% CI [%.3e, %.3e]\n", name, q[2], q[1], q[3])
    p1 = plot(obj.Output[:, i], title="chain: $name", label=""); vline!(p1,[burn_in],lc=:red,ls=:dash,label="")
    p2 = histogram(c, title="posterior: $name", label=""); vline!(p2,[QoI_dict[name][2]],lc=:green,lw=3,label="prior center")
    savefig(p1, joinpath(figdir, "venting_$(name)_chain.png"))
    savefig(p2, joinpath(figdir, "venting_$(name)_posterior.png"))
end

# posterior-median fit overlay (daytime)
_, T2m, RHm, _ = forward_map(med, QoI, base, days, tspan, teval; u0=copy(u0))
ph = (predtimes .- t0)./3600.0
pT = plot(ph, Tmeas .- 273.15, lc=:black, lw=3, label="measured T", xlabel="hours", ylabel="T air (C)")
plot!(pT, ph, collect(T2m) .- 273.15, lc=:orange, lw=2, label="model (posterior median)")
vspan!(pT, [DAY_LO, DAY_HI], alpha=0.08, lc=:blue, label="scored window")
savefig(pT, joinpath(figdir, "venting_fit_T.png"))
println("\nFigures: venting_energy.png, venting_*_posterior.png, venting_fit_T.png in ", figdir)
println("""
Read: if nu4's posterior mass sits near ~1e-4 and well below the old 1e-2, the data
confirm the physical leakage under the corrected physics. A broad nu4 posterior that
simply excludes large values is still a valid result (on an open-vent day the roof
vent does the cooling, so the data only need leakage to be small). gamma4 should be
well identified by RH.
""")
