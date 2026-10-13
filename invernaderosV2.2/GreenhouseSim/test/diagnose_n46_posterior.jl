# =============================================================================
# diagnose_n46_posterior.jl -- what the n46 posterior does, and what it costs to keep it physical (Note 49)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     julia --project=. test/diagnose_n46_posterior.jl
#
# Uses the n46 chains (final file if present, else the checkpoint), burn-in 30 000, pooled over the seeds found.
# Runs (full season, closed loop, n46 structure):
#   P     pooled posterior mean
#   P1    P with f_mass = 0.15 (cap: inside solar absorption ~0.8 of outside)
#   P2    P1 with k_deep = 2 (no large deep store)
#   P0    P with f_mass = 0 (what the sunlit mass contributes)
# For each: hourly RMSE, night bias by period, day bias by hour of day, VP, vent bias, and a "sensor test":
#   regression of the daytime residual (measured - model) on outside radiation, resid = a + b * Iglob.
#   A clear positive b (K per 100 W m-2) with a small a means the error scales with sunlight: either missing
#   solar heat in the model or a sun-heated sensor. The two cannot be told apart from air T alone.
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n46"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

const BURN46 = 30000
parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n46_seed$(s).jld2")
    fc = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n46_seed$(s)_ckpt.jld2")
    src = isfile(f) ? f : isfile(fc) ? fc : ""
    isempty(src) && continue
    ch = JLD2.load(src)["chain"]
    size(ch, 1) > BURN46 + 1000 || continue
    push!(parts, ch[BURN46+1:end, 1:NM])
    @printf("seed %d: %s, %d post-burn steps\n", s, basename(src), size(ch, 1) - BURN46)
end
isempty(parts) && error("no n46 chains found in $POSTDIR")
pm = vec(mean(reduce(vcat, parts), dims = 1))
println("\nPooled posterior mean:")
for (q, v) in zip(CQ_MODEL, pm); @printf("  %-14s %10.4g\n", q, v); end
ix(q) = findfirst(==(q), CQ_MODEL)
@printf("\nImplied deep-soil capacity k_deep * sd_tau = %.2g J m-2 K-1 (~%.1f m of soil at 2e6 J m-3 K-1)\n",
        pm[ix("k_deep")] * params_cl.sd_tau, pm[ix("k_deep")] * params_cl.sd_tau / 2e6)
@printf("Mass time constant C_mass / h_mass = %.0f min\n", pm[ix("C_mass")] / pm[ix("h_mass")] / 60)

function setx(x, pairs...)
    y = copy(x); for (q, v) in pairs; y[ix(q)] = v; end; y
end
runs = [("P  posterior mean", pm),
        ("P1 f_mass 0.15", setx(pm, "f_mass" => 0.15)),
        ("P2 P1 + k_deep 2", setx(pm, "f_mass" => 0.15, "k_deep" => 2.0)),
        ("P0 f_mass 0", setx(pm, "f_mass" => 0.0))]

okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g); okU = .!isnan.(U8obs); okU[1] = false
inset(ds) = in.(cd_, Ref(Set(ds)))
ALLD = inset(collect(0:114)); COM = inset(COMMON_TEST)
VP(T, RH) = RH ./ 100 .* Pws.(T) ./ 100; VPo = VP(Tobs_g, RHobs_g)
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
function hourly_rmse(v, o, mask)
    e = Float64[]
    for cols in values(HB); k = [c for c in cols if mask[c]]; length(k) >= 6 && push!(e, mean(v[k]) - mean(o[k])); end
    isempty(e) ? NaN : sqrt(mean(abs2, e))
end
PER = [(0, 8), (9, 48), (49, 100), (101, 109), (110, 114)]
hrc = [floor(Int, mod(((c - 2) % NPD) * 300 / 3600, 24.0)) for c in 1:NCOL]; hrc[1] = 0
HB_ = ((6, 9), (9, 12), (12, 15), (15, 18), (18, 21))
Igc = A(Ig_v).(predtimes)

println("\n", "="^170)
@printf("  %-20s | %6s %6s | %s | %s | %6s %6s %6s %7s | %s\n", "run", "hRMS", "common",
        join([@sprintf("%6s", "n$(a)") for (a, b) in PER], " "), join([@sprintf("%6s", "d$(a)-$(b)") for (a, b) in HB_], " "),
        "T day", "VPday", "VPnit", "U8 day", "sensor test: resid = a + b*I  (a K, b K/100 W m-2)")
println("  ", "-"^168)
for (nm_, x) in runs
    T, RH, U = full_eval(x)
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    m = okT .& dayt; r = Tobs_g[m] .- T[m]; I = Igc[m]
    X = hcat(ones(length(I)), I ./ 100); ab = X \ r
    @printf("  %-20s | %6.3f %6.3f | %s | %s | %+6.2f %+6.2f %+6.2f %+7.3f | a %+.2f  b %+.3f\n", nm_,
            hourly_rmse(T, Tobs_g, okT .& ALLD), hourly_rmse(T, Tobs_g, okT .& COM),
            join([@sprintf("%+6.2f", v) for v in nb], " "), join([@sprintf("%+6.2f", v) for v in hb], " "),
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& nite),
            mb(U, U8obs, okU .& ALLD .& dayt), ab[1], ab[2])
end
println("""
  How to read it:
   - P vs P1/P2: the cost, in RMSE and bias, of keeping the solar input and the deep store physical.
   - P0: how much the sunlit mass does at the posterior.
   - Sensor test: b > 0 means the measured air is warmer than the model in proportion to sunlight. With
     f_mass already at its upper bound in P, a remaining b > 0 points to either more solar heat than physically
     plausible or a sun-heated sensor; the AGC-2018 sensor documentation would decide.""")
