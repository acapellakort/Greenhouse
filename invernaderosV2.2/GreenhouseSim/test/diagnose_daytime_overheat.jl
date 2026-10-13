# =============================================================================
# diagnose_daytime_overheat.jl -- why does the model overheat on closed-vent sunny days?
# =============================================================================
# Run from the GreenhouseSim folder (~10 s of simulation, no MCMC, no src/ changes):
#     julia --project=. test/diagnose_daytime_overheat.jl
#
# Spike-day inspection (Note 32-d) showed: no missing vent data, no blackout screen by day;
# the real greenhouse holds a flat 26-27 C at 700-870 W/m2 with vents 5-25 %, while the
# model swings 10-20 K with every small change in the logged opening. Candidates, each
# changed ALONE from N29 with U8_min_leak = 0 (all are input/parameter changes, no engine edit):
#   D0 base        N29, U8_min_leak = 0
#   D1 vent^0.5    effective opening = U8^0.5 (small openings ventilate more than linear;
#                  real vent flow vs opening % is concave). Night vents are mostly ~0, so
#                  this mostly acts by day -- unlike the global U8_min_leak.
#   D2 vent^0.3    stronger concave map
#   D3 transp max  gamma4 = 20 (prior lower bound) + derived observed LAI(t): maximum latent cooling
#   D4 tau1 x0.85  15 % less solar transmission (dirty/whitewashed cover)
#   D5 leak 0.05   the Note 29 setting, for reference
# Metrics by window: days 0-10 (vents open), 11-41 (spike period), 42-114 (holdout).
# "Roughness" = RMS of 5-min changes by day: how jumpy the model is vs the measurement.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = get(ENV, "N_FULL", "115")
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const DO_FIG = joinpath(FIGDIR, "AiCU_daytime_overheat")

lai_d = vcat(0.0, Float64.(cm_days_all)); lai_v = vcat(0.5, cm_L_all)
lai_of_day(x) = x <= lai_d[1] ? lai_v[1] : x >= lai_d[end] ? lai_v[end] :
    (i = searchsortedlast(lai_d, x); lai_v[i] + (x - lai_d[i]) / (lai_d[i+1] - lai_d[i]) * (lai_v[i+1] - lai_v[i]))
weather_lai = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, _t -> 18.0 + 273.15,
                            A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), tcal -> lai_of_day((tcal - t0) / 86400.0))
ctrl_u8(f) = ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(f.(U8v)), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipev))
base_v05 = SimBase(params, weather, ctrl_u8(u -> sqrt(u)), t0)
base_v03 = SimBase(params, weather, ctrl_u8(u -> u^0.3), t0)
base_lai = SimBase(params, weather_lai, controls, t0)

nV = QoI_V; θ0 = [4.46e-5, 44.0, 24.5, 0.0]
runs = [("D0 base",       nV,               θ0,                           base),
        ("D1 vent^0.5",   nV,               θ0,                           base_v05),
        ("D2 vent^0.3",   nV,               θ0,                           base_v03),
        ("D3 transp max", nV,               [θ0[1], 20.0, θ0[3], 0.0],    base_lai),
        ("D4 tau1 x0.85", vcat(nV, "tau1"), vcat(θ0, 0.85 * params.tau1), base),
        ("D5 leak 0.05",  nV,               [θ0[1:3]; 0.05],              base)]
T2s = Vector{Vector{Float64}}(); RHs = Vector{Vector{Float64}}()
println("Running $(length(runs)) full-season forward simulations ...")
for (name, nm, x, b) in runs
    _, T2, RH, _ = forward_map(x, nm, b, days, tspan, teval; u0 = copy(u0))
    length(T2) == NCOL || error("$name: solver stopped early")
    push!(T2s, collect(T2)); push!(RHs, collect(RH))
end

coldays = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
okT  = .!isnan.(Tobs_g); okT[1] = false
okR  = .!isnan.(RHobs_g); okR[1] = false
dayt = A(Ig_v).(predtimes) .> 50.0
nite = A(Ig_v).(predtimes) .< 1.0
win(a, b) = (coldays .>= a) .& (coldays .<= b)
W = [("0-10", win(0, 10)), ("11-41", win(11, 41)), ("42-114", win(42, N_DAYS - 1))]
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
rmsem(v, o, m) = any(m) ? sqrt(mean((v[m] .- o[m]) .^ 2)) : NaN

# daily max overshoot (model daytime max - measured daytime max) and roughness
function daystats(T, d)
    m = okT .& dayt .& (coldays .== d)
    count(m) < 12 && return (NaN, NaN, NaN)
    ix = findall(m)
    rough(v) = sqrt(mean(abs2, diff(v[ix])))
    return (maximum(T[m]) - maximum(Tobs_g[m]), rough(T), rough(Tobs_g))
end

println("\n", "="^118)
println("  DAYTIME (Iglob > 50) temperature [K]               | spike period 11-41               | DAY RH bias [%] | NIGHT bias [K]")
@printf("  %-14s %8s %8s %8s %8s | %7s %7s %8s %7s | %7s %7s | %7s %7s\n", "run", "bias", "bias", "RMSE", "bias",
        "mean", "days", "rough", "rough", "0-41", "42-114", "0-41", "42-114")
@printf("  %-14s %8s %8s %8s %8s | %7s %7s %8s %7s | %7s %7s | %7s %7s\n", "", "0-10", "11-41", "11-41", "42-114",
        "ovrsht", ">8K", "model", "meas", "", "", "", "")
println("  ", "-"^116)
for (k, (name, _, _, _)) in enumerate(runs)
    T = T2s[k]; R = RHs[k]
    ds = [daystats(T, d) for d in 11:41]; ov = filter(!isnan, first.(ds))
    @printf("  %-14s %+8.2f %+8.2f %8.2f %+8.2f | %+7.2f %7d %8.2f %7.2f | %+7.1f %+7.1f | %+7.2f %+7.2f\n", name,
            mb(T, Tobs_g, okT .& dayt .& W[1][2]), mb(T, Tobs_g, okT .& dayt .& W[2][2]), rmsem(T, Tobs_g, okT .& dayt .& W[2][2]),
            mb(T, Tobs_g, okT .& dayt .& W[3][2]), mean(ov), count(>(8), ov),
            mean(filter(!isnan, getindex.(ds, 2))), mean(filter(!isnan, getindex.(ds, 3))),
            mb(R, RHobs_g, okR .& dayt .& win(0, 41)), mb(R, RHobs_g, okR .& dayt .& W[3][2]),
            mb(T, Tobs_g, okT .& nite .& win(0, 41)), mb(T, Tobs_g, okT .& nite .& W[3][2]))
end
println("""
  Read: the right fix lowers the 11-41 overshoot and the model roughness toward the measured roughness
  WITHOUT pushing days 0-10 or the nights colder. D3 (more transpiration) should also raise RH; D4 (less
  sun) should lower it; D1/D2 (more ventilation at small openings) should lower RH and act mainly by day.""")

# figure: days 24, 30, 41 at 5-min, all runs
SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "21,24,30,41"), ","))
ps = Any[]
for d in SHOWD
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    p = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 3, label = "measured", title = "day $d", xlabel = "hour",
             ylabel = "T air (°C)", legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6)
    for (k, (name, _, _, _)) in enumerate(runs)
        plot!(p, h, T2s[k][c] .- 273.15, lw = 1.2, label = name)
    end
    push!(ps, p)
end
savefig(plot(ps..., layout = (1, length(ps)), size = (420 * length(ps), 380)), DO_FIG * ".png")
println("Figure: ", DO_FIG, ".png")
