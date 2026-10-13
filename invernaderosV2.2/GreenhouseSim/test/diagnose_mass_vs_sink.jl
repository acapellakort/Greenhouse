# =============================================================================
# diagnose_mass_vs_sink.jl -- thermal MASS (capacitance) or latent SINK (transpiration)?
# =============================================================================
# Run from the GreenhouseSim folder (~10 s of simulation, no MCMC, no src/ changes):
#     julia --project=. test/diagnose_mass_vs_sink.jl
#
# One change at a time from N29 with U8_min_leak = 0 (existing parameters/inputs only):
#   M0 base
#   M1 canopy mass x5     alpha1 x5   (canopy capacity 3 -> 15 kJ/m2/K at LAI 2)
#   M2 canopy mass x20    alpha1 x20  (-> 60 kJ/m2/K; canopy-air time constant ~50 min)
#   M3 sink: gamma4 = 20  minimum stomatal resistance at the prior bound, LAI = 2.0
#   M4 LAI(t) only        derived observed LAI (0.5 -> 2.5), gamma4 = 44
#   M5 mass x20 + sink    M2 and M3 together (are the effects additive?)
#   F  sensor filter      M0 output passed through a first-order 10-min filter (post hoc;
#                         an aspirated sensor box lags the air -- can it explain the roughness?)
# Expected signatures:
#   mass  -> roughness and overshoot down, daily-mean daytime bias ~unchanged, VP ~unchanged
#   sink  -> daytime mean bias down AND daytime vapour pressure up (toward measured)
# Vapour pressure (VP) is used instead of RH because RH is confounded with temperature.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = get(ENV, "N_FULL", "115")
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const MS_FIG  = joinpath(FIGDIR, "AiCU_mass_vs_sink")
const TAU_SEN = parse(Float64, get(ENV, "TAU_SENSOR_MIN", "10")) * 60

lai_d = vcat(0.0, Float64.(cm_days_all)); lai_v = vcat(0.5, cm_L_all)
lai_of_day(x) = x <= lai_d[1] ? lai_v[1] : x >= lai_d[end] ? lai_v[end] :
    (i = searchsortedlast(lai_d, x); lai_v[i] + (x - lai_d[i]) / (lai_d[i+1] - lai_d[i]) * (lai_v[i+1] - lai_v[i]))
base_lai = SimBase(params, WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, _t -> 18.0 + 273.15,
                   A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), tcal -> lai_of_day((tcal - t0) / 86400.0)), controls, t0)

nV = QoI_V; θ0 = [4.46e-5, 44.0, 24.5, 0.0]; a1 = params.alpha1
runs = [("M0 base",          nV,                 θ0,                                     base),
        ("M1 mass x5",       vcat(nV, "alpha1"), vcat(θ0, 5a1),                          base),
        ("M2 mass x20",      vcat(nV, "alpha1"), vcat(θ0, 20a1),                         base),
        ("M3 sink g4=20",    nV,                 [θ0[1], 20.0, θ0[3], 0.0],              base),
        ("M4 LAI(t)",        nV,                 θ0,                                     base_lai),
        ("M5 mass20+sink",   vcat(nV, "alpha1"), [θ0[1], 20.0, θ0[3], 0.0, 20a1],        base)]
T2s = Vector{Vector{Float64}}(); RHs = Vector{Vector{Float64}}()
println("Running $(length(runs)) full-season forward simulations ...")
for (name, nm, x, b) in runs
    t = @elapsed (out = forward_map(x, nm, b, days, tspan, teval; u0 = copy(u0)))
    length(out[2]) == NCOL || error("$name: solver stopped early")
    push!(T2s, collect(out[2])); push!(RHs, collect(out[3]))
    @printf("  %-16s %.1f s\n", name, t)
end
# post-hoc sensor filter on M0 (exact first-order discretisation on the 5-min grid)
af = exp(-300.0 / TAU_SEN)
filt(x) = (y = similar(x); y[1] = x[1]; for k in 2:length(x); y[k] = af * y[k-1] + (1 - af) * x[k]; end; y)
push!(runs, (@sprintf("F  M0+sensor %.0fmin", TAU_SEN / 60), nV, θ0, base))
push!(T2s, filt(T2s[1])); push!(RHs, RHs[1])

coldays = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
okT  = .!isnan.(Tobs_g); okT[1] = false
okV  = okT .& .!isnan.(RHobs_g)
dayt = A(Ig_v).(predtimes) .> 50.0
nite = A(Ig_v).(predtimes) .< 1.0
win(a, b) = (coldays .>= a) .& (coldays .<= b)
VPobs = RHobs_g ./ 100 .* Pws.(Tobs_g) ./ 100           # hPa
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN

function spike_stats(T)
    ov = Float64[]; rmod = Float64[]; rmea = Float64[]
    for d in 11:41
        m = okT .& dayt .& (coldays .== d); count(m) < 12 && continue
        ix = findall(m)
        push!(ov, maximum(T[m]) - maximum(Tobs_g[m]))
        push!(rmod, sqrt(mean(abs2, diff(T[ix])))); push!(rmea, sqrt(mean(abs2, diff(Tobs_g[ix]))))
    end
    return mean(ov), count(>(8), ov), mean(rmod), mean(rmea)
end

println("\n", "="^124)
println("  Temperature [K]: DAYTIME (Iglob>50) bias              | spike period 11-41               | VAPOUR PRESSURE bias [hPa] | NIGHT T bias")
@printf("  %-20s %7s %7s %7s | %7s %5s %7s %6s | %7s %7s %7s | %7s %7s\n", "run", "0-10", "11-41", "42-114",
        "ovrsht", ">8K", "rough", "meas", "day0-41", "day42+", "night", "0-41", "42-114")
println("  ", "-"^122)
for (k, r) in enumerate(runs)
    T = T2s[k]; VP = RHs[k] ./ 100 .* Pws.(T2s[k]) ./ 100
    ov, n8, rm_, rmeas = spike_stats(T)
    @printf("  %-20s %+7.2f %+7.2f %+7.2f | %+7.2f %5d %7.2f %6.2f | %+7.2f %+7.2f %+7.2f | %+7.2f %+7.2f\n", r[1],
            mb(T, Tobs_g, okT .& dayt .& win(0, 10)), mb(T, Tobs_g, okT .& dayt .& win(11, 41)),
            mb(T, Tobs_g, okT .& dayt .& win(42, N_DAYS - 1)), ov, n8, rm_, rmeas,
            mb(VP, VPobs, okV .& dayt .& win(0, 41)), mb(VP, VPobs, okV .& dayt .& win(42, N_DAYS - 1)),
            mb(VP, VPobs, okV .& nite), mb(T, Tobs_g, okT .& nite .& win(0, 41)), mb(T, Tobs_g, okT .& nite .& win(42, N_DAYS - 1)))
end
println("""
  F row: temperature only (vapour pressure columns repeat M0, the filter is not applied to RH).
  Mass signature: overshoot and roughness fall, daytime MEAN bias and VP barely move.
  Sink signature: daytime mean bias falls AND daytime VP rises toward the measurement.""")

SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "21,24,30,41"), ","))
ps = Any[]
for d in SHOWD
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    p = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 3, label = "measured", title = "day $d", xlabel = "hour",
             ylabel = "T air (°C)", legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6)
    for (k, r) in enumerate(runs)
        k in (1, 2, 3, 4, 6, 7) && plot!(p, h, T2s[k][c] .- 273.15, lw = k == 1 ? 1 : 1.6, label = r[1])
    end
    push!(ps, p)
end
savefig(plot(ps..., layout = (1, length(ps)), size = (420 * length(ps), 380)), MS_FIG * ".png")
println("Figure: ", MS_FIG, ".png")
