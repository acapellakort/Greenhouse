# =============================================================================
# diagnose_n46_screen.jl -- vents, pipe heat (AGC energy formula) and physical caps around the n46 posterior (Note 49)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     julia --project=. test/diagnose_n46_screen.jl
#
# Uses the n46 chains (final file if present, else the checkpoint), burn-in 30 000, pooled over the seeds found.
# Runs (full season, closed loop, n46 structure), starting from the pooled posterior mean P:
#   V1/V2  vent discharge nu1 x0.8 / x0.6        (day cold AND dry at P -> too much exchange per opening?)
#   W1     wind pressure coefficient nu2 x0.5    (wind term dominates the roof-vent flow)
#   H1     pipe emission matched to the AGC heating-energy formula: Heating_Energy = 3.37 * (PipeLow - Tair)^1.32
#          W m-2 (R2 0.993, data-only check); model convection 0.638*gamma1 + ~0.2 radiation -> gamma1 = 4.9
#   H2     H1 + k_deep 2 (no deep store)          H3  H2 + f_mass 0.15 (physical solar share)
#   H4     H3 + nu1 x0.7                          -- the "all physical" candidate
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
nu1P = pm[ix("nu1")]
runs = [("P  posterior mean", pm),
        ("V1 nu1 x0.8", setx(pm, "nu1" => 0.8nu1P)),
        ("V2 nu1 x0.6", setx(pm, "nu1" => 0.6nu1P)),
        ("H1 gamma1 4.9 (AGC)", setx(pm, "gamma1" => 4.9)),
        ("H2 H1 + k_deep 2", setx(pm, "gamma1" => 4.9, "k_deep" => 2.0)),
        ("H3 H2 + f_mass 0.15", setx(pm, "gamma1" => 4.9, "k_deep" => 2.0, "f_mass" => 0.15)),
        ("H4 H3 + nu1 x0.7", setx(pm, "gamma1" => 4.9, "k_deep" => 2.0, "f_mass" => 0.15, "nu1" => 0.7nu1P))]
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
# W1: wind pressure coefficient is not a calibrated parameter -> run with a modified params_cl
let pW = update_params(update_params(params_cl, CQ_MODEL, pm), ["nu2"], [0.5 * params_cl.nu2])
    us = simdays(pW, days, u0_cl, RELTOL); T = us[2, :]; RH = rhf.(us[2, :], us[4, :]); U = clamp.(us[7, :], 0.0, 1.0)
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    m = okT .& dayt; r = Tobs_g[m] .- T[m]; I = Igc[m]; ab = hcat(ones(length(I)), I ./ 100) \ r
    @printf("  %-20s | %6.3f %6.3f | %s | %s | %+6.2f %+6.2f %+6.2f %+7.3f | a %+.2f  b %+.3f\n", "W1 nu2 x0.5 (wind)",
            hourly_rmse(T, Tobs_g, okT .& ALLD), hourly_rmse(T, Tobs_g, okT .& COM),
            join([@sprintf("%+6.2f", v) for v in nb], " "), join([@sprintf("%+6.2f", v) for v in hb], " "),
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& nite),
            mb(U, U8obs, okU .& ALLD .& dayt), ab[1], ab[2])
end
println("""
  How to read it:
   - V/W rows: if less exchange per opening warms AND humidifies the day while U8 day stays near 0, the vent
     flow per opening was too large (compensated before by gamma4 and the solar share).
   - H rows: with the pipe heat the AGC accounting implies, can the deep store (k_deep) and the oversized solar
     share be dropped without cold nights? H4 is the all-physical candidate for a re-calibration.
   - Sensor test: b > 0 means the measured air is warmer than the model in proportion to sunlight. With
     f_mass already at its upper bound in P, a remaining b > 0 points to either more solar heat than physically
     plausible or a sun-heated sensor; the AGC-2018 sensor documentation would decide.""")
