# =============================================================================
# diagnose_posterior.jl -- season breakdown of a calibrated posterior (n46, n49a, n49b) (Notes 49, 52)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     PRIORS=n49b julia --project=. test/diagnose_posterior.jl      (PRIORS = n46 | n49a | n49b)
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

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n46")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

const BURN46 = 30000
parts = Vector{Matrix{Float64}}()
const SEEDS_ = parse.(Int, split(get(ENV, "SEEDS", "1,2"), ","))   # Note 59: SEEDS=2 to evaluate one mode of a bimodal posterior
for s in SEEDS_
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2")
    fc = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s)_ckpt.jld2")
    src = isfile(f) ? f : isfile(fc) ? fc : ""
    isempty(src) && continue
    ch = JLD2.load(src)["chain"]
    size(ch, 1) > BURN46 + 1000 || continue
    push!(parts, ch[BURN46+1:end, 1:NM])
    @printf("seed %d: %s, %d post-burn steps\n", s, basename(src), size(ch, 1) - BURN46)
end
isempty(parts) && error("no $(PRIORS) chains found in $POSTDIR")
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
runs = Any[("P  posterior mean", pm)]
if "dday" in CQ_MODEL
    push!(runs, ("P  air T only (dday 0)", setx(pm, "dday" => 0.0)))     # what the physical air node does
else
    push!(runs, ("P1 f_mass 0.15", setx(pm, "f_mass" => 0.15)))
end
if "s_rad" in CQ_MODEL     # n61: what the physical air node does without the sun-proportional measurement term
    push!(runs, ("P  air T only (s_rad 0)", setx(pm, "s_rad" => 0.0)))
    push!(runs, ("P  air cap x1, s_rad 0", setx(pm, "s_rad" => 0.0, "air_cap_mult" => 1.0)))
end
# Note 60 (SOLARGRID=1): forward screen of the afternoon solar deficit at the posterior mean
#   M*   measurement model: compared T = air T + s * Iglob/100 (sensor heated by sun; post hoc, no re-simulation)
#   F*   f_mass above its 0.2 cap (more absorbed sun reaching the air via the internal mass)
#   G*   gamma4 higher (less daytime transpiration / latent cooling)
#   FG   both
const SOLARGRID = get(ENV, "SOLARGRID", "0") == "1"
if SOLARGRID
    runs = Any[("P  posterior mean", pm)]
    for s_ in (0.15, 0.3); push!(runs, ("M  sensor s=$(s_)", pm)); end
    for f in (0.25, 0.3, 0.4); push!(runs, ("F  f_mass $(f)", setx(pm, "f_mass" => f))); end
    for g in (110.0, 150.0); push!(runs, ("G  gamma4 $(Int(g))", setx(pm, "gamma4" => g))); end
    push!(runs, ("FG f_mass 0.3 g4 110", setx(pm, "f_mass" => 0.3, "gamma4" => 110.0)))
end
# Note 60 (LAGGRID=1): first-order SENSOR LAG on the compared air T (post hoc filter of the 5-min model series,
#   tau in minutes; no re-simulation). A lagged sensor makes the measured diurnal cycle trail the air and damps
#   within-hour swings -- the two signatures of Notes 58/59 (morning warm / afternoon cold; within-hour SD x2).
const LAGGRID = get(ENV, "LAGGRID", "0") == "1"
if LAGGRID
    runs = Any[("P  posterior mean", pm)]
    for tau in (15, 30, 60, 90, 120); push!(runs, ("L  lag $(tau) min", pm)); end
    push!(runs, ("LM lag 60 + s=0.15", pm))
    push!(runs, ("LM lag 30 + s=0.15", pm))
    for m in (2.0, 4.0, 8.0); push!(runs, ("A  air cap x$(Int(m))", pm)); end   # Note 61: physical alternative to a sensor lag
end
function lagT(T, tau_min)
    a = 300.0 / (60.0 * tau_min); y = copy(T)
    for k in 3:length(T); y[k] = y[k-1] + a * (T[k] - y[k-1]); end   # column 1 is the initial state
    y
end
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
TP = Float64[]
Igc_ = A(Ig_v).(predtimes)
for (nm_, x) in runs
    if startswith(nm_, "A  air cap x")       # re-simulate with a larger effective air heat capacity
        m_ = parse(Float64, split(nm_, "x")[end])
        us_ = simdays(update_params(update_params(params_cl, CQ_MODEL, x), ["air_cap_mult"], [m_]), days, u0_cl, RELTOL)
        T, RH = us_[2, :], rhf.(us_[2, :], us_[4, :]); U = AU.(predtimes)
    else
        T, RH, U = full_eval(x)
    end
    Tair = copy(T)          # VP always from the air node
    if startswith(nm_, "M  sensor s=")          # measurement model: compare T + s*I/100 (RH kept at the air node)
        s_ = parse(Float64, split(nm_, "=")[2]); T = T .+ s_ .* Igc_ ./ 100
    end
    if startswith(nm_, "L  lag ")
        T = lagT(T, parse(Float64, split(nm_)[3]))
    elseif startswith(nm_, "LM ")
        T = lagT(T, parse(Float64, split(nm_)[3])) .+ 0.15 .* Igc_ ./ 100
    end
    if LAGGRID     # within-hour SD of T (model as compared / observed), common test, day and night
        wsd(v, m) = (e = Float64[]; for cols in values(HB); k = [c for c in cols if m[c]]; length(k) >= 6 && push!(e, std(v[k])); end; mean(e))
        @printf("  %-20s   within-hour SD  day %.3f / %.3f   night %.3f / %.3f\n", nm_, wsd(T, okT .& COM .& dayt),
                wsd(Tobs_g, okT .& COM .& dayt), wsd(T, okT .& COM .& nite), wsd(Tobs_g, okT .& COM .& nite))
    end
    isempty(TP) && append!(TP, T)
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    m = okT .& dayt; r = Tobs_g[m] .- T[m]; I = Igc[m]
    X = hcat(ones(length(I)), I ./ 100); ab = X \ r
    @printf("  %-20s | %6.3f %6.3f | %s | %s | %+6.2f %+6.2f %+6.2f %+7.3f | a %+.2f  b %+.3f\n", nm_,
            hourly_rmse(T, Tobs_g, okT .& ALLD), hourly_rmse(T, Tobs_g, okT .& COM),
            join([@sprintf("%+6.2f", v) for v in nb], " "), join([@sprintf("%+6.2f", v) for v in hb], " "),
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(VP(Tair, RH), VPo, okV .& ALLD .& dayt), mb(VP(Tair, RH), VPo, okV .& ALLD .& nite),
            mb(U, U8obs, okU .& ALLD .& dayt), ab[1], ab[2])
end
# Note 55: Note 45b criteria per season period, at the posterior mean (first run)
println("\n  Note 45b criteria by period (posterior mean; model - measured, K):")
@printf("  %-10s %8s %8s %8s   %s\n", "period", "24-h", "day", "night", "crop: |24-h| <= 0.5 ; physics: |day|,|night| <= 0.75")
for (a, b) in PER
    md = inset(collect(a:b))
    v24 = mb(TP, Tobs_g, okT .& md); vd = mb(TP, Tobs_g, okT .& md .& dayt); vn = mb(TP, Tobs_g, okT .& md .& nite)
    @printf("  d%-3d-%-4d %+8.2f %+8.2f %+8.2f   %s %s\n", a, b, v24, vd, vn,
            abs(v24) <= 0.5 ? "crop ok  " : "crop FAIL", (abs(vd) <= 0.75 && abs(vn) <= 0.75) ? "physics ok" : "physics FAIL")
end
println("""
  How to read it:
   - P vs P1/P2: the cost, in RMSE and bias, of keeping the solar input and the deep store physical.
   - P0: how much the sunlit mass does at the posterior.
   - Sensor test: b > 0 means the measured air is warmer than the model in proportion to sunlight. With
     f_mass already at its upper bound in P, a remaining b > 0 points to either more solar heat than physically
     plausible or a sun-heated sensor; the AGC-2018 sensor documentation would decide.""")
