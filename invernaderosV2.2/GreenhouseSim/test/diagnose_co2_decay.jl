# =============================================================================
# diagnose_co2_decay.jl -- night CO2 decay test of the air exchange, independent of the daytime level (Note 65)
# =============================================================================
# Run from the GreenhouseSim folder (~1 min per model, forward only, no MCMC):
#     PRIORS=n54 julia --project=. test/diagnose_co2_decay.jl      (closed loop: model vents)
#     PRIORS=n61 julia --project=. test/diagnose_co2_decay.jl      (open loop: logged vents)
#
# Why: in Note 64 the night CO2 level was inherited from the daytime dosing deficit, so it said nothing about
# ventilation. Here the level is removed: on night hours with lamps OFF and no dosing (U10 = 0) both the measured
# and the simulated CO2 are fitted, night by night, to
#         dC/dt = R - k (C - C_out)            (hourly means; k = air exchange [1/h], R = net source)
# k and R come from the SHAPE of each night's curve, not its level -- equivalent to resetting the model to the
# measured CO2 at dusk, without a re-simulation. C_out = 834.7 mg m-3, the model's constant outdoor value.
#   1. by season period: pooled within-night k (night-specific intercepts) for measured vs model, and the median
#      per-night R (ppm/h) at that k.
#   2. nights with vents CLOSED (logged opening < 0.02: leakage only) vs OPEN, and k vs wind on closed nights.
#   3. night-by-night agreement: correlation of per-night k (model vs measured).
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n61")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2"); fc = replace(f, ".jld2" => "_ckpt.jld2")
    src = isfile(f) ? f : isfile(fc) ? fc : ""; isempty(src) && continue
    ch = JLD2.load(src)["chain"]; size(ch, 1) > 31000 || continue
    push!(parts, ch[30001:end, 1:NM]); println("seed $s: ", basename(src))
end
isempty(parts) && error("no $(PRIORS) chains in $POSTDIR")
pmean = vec(mean(reduce(vcat, parts), dims = 1))
usd = simdays(update_params(params_cl, CQ_MODEL, pmean), days, u0_cl, RELTOL); usd === nothing && error("solve failed")

C_OUTD = 834.7                                   # model outdoor CO2 (mg m-3)
toppm  = 293.15 / 536.4                          # mg m-3 -> ppm at 20 C (for R only; k is unit-free)
Cmod = usd[3, :]; Cobs = AC.(predtimes)          # both mg m-3 (measured converted with measured T at load)
lampon = A(U12v).(predtimes) .> 0.5; dosing = A(U10v).(predtimes) .> 0.0
Vlog = AU.(predtimes); Wnd = A(Wind_v).(predtimes)
Vmod = OLP ? Vlog : clamp.(usd[7, :], 0.0, 1.0)
use = .!isnan.(Cobs) .& nite .& .!lampon .& .!dosing; use[1] = false

hix = [floor(Int, (predtimes[c] - t0) / 3600) for c in 1:NCOL]
HH = Dict{Int,Vector{Int}}()
for c in 2:NCOL; use[c] && push!(get!(HH, hix[c], Int[]), c); end
nhrs = sort([h for (h, v) in HH if length(v) >= 10])
nights = Vector{Vector{Int}}(); cur = Int[]
for h in nhrs
    global cur
    if !isempty(cur) && h != cur[end] + 1; push!(nights, cur); cur = Int[]; end
    push!(cur, h)
end
isempty(cur) || push!(nights, cur)
nights = filter(r -> length(r) >= 5, nights)

hmean(v, h) = mean(v[HH[h]])
function xy(C, r)                 # x = mid-hour (C - C_out), y = hourly change (mg m-3 per h)
    c = [hmean(C, h) for h in r]
    return (c[1:end-1] .+ c[2:end]) ./ 2 .- C_OUTD, diff(c)
end
function kfit(C, rs)              # pooled within-night k: night-specific intercepts (demeaned per night)
    sxy = 0.0; sxx = 0.0
    for r in rs; x, y = xy(C, r); x = x .- mean(x); y = y .- mean(y); sxy += sum(x .* y); sxx += sum(x .^ 2); end
    return sxx > 0 ? -sxy / sxx : NaN
end
Rnight(C, r, k) = (xy_ = xy(C, r); mean(xy_[2]) + k * mean(xy_[1]))
function k1(C, r)                 # per-night least squares (k, R)
    x, y = xy(C, r); X = hcat(ones(length(x)), -x); b = X \ y; return b[2], b[1]
end
nday(r)  = max(0, floor(Int, (r[1] - 12) / 24))       # day of the evening on which the night starts
nvent(r) = mean(vcat([Vlog[HH[h]] for h in r]...))
nwind(r) = mean(vcat([Wnd[HH[h]] for h in r]...))
nrange(C, r) = (c = [hmean(C, h) for h in r]; maximum(c) - minimum(c))

println("\n", "="^112)
@printf("  %s  night CO2 decay test (lamps off, no dosing), posterior mean, %s\n", PRIORS, OLP ? "OPEN loop (logged vents)" : "CLOSED loop (model vents)")
@printf("  %d nights with >= 5 consecutive usable hours (%d hours)\n", length(nights), sum(length, nights))
println("\n  1. Air exchange k (1/h) from the night decay shape; R = net night source at that k (ppm/h)")
@printf("  %-10s | %6s | %8s %8s %7s | %8s %8s | %8s %8s | %6s %6s\n", "period", "nights", "k obs", "k model", "ratio",
        "R obs", "R model", "C obs", "C model", "vent", "wind")
PERD = [(0, 8), (9, 48), (49, 100), (101, 109), (110, 114), (0, 114)]
for (a, b) in PERD
    rs = [r for r in nights if a <= nday(r) <= b]; isempty(rs) && continue
    ko = kfit(Cobs, rs); km = kfit(Cmod, rs)
    @printf("  d%-3d-%-4d | %6d | %8.3f %8.3f %7.2f | %+8.1f %+8.1f | %8.0f %8.0f | %6.3f %6.2f\n", a, b, length(rs), ko, km, km / ko,
            median([Rnight(Cobs, r, ko) for r in rs]) * toppm, median([Rnight(Cmod, r, km) for r in rs]) * toppm,
            mean([hmean(Cobs, h) for r in rs for h in r]) * toppm, mean([hmean(Cmod, h) for r in rs for h in r]) * toppm,
            mean(nvent.(rs)), mean(nwind.(rs)))
end
println("\n  2. Vent state (logged opening averaged over the night)")
for (lab, f) in (("closed  (< 0.02, leakage only)", r -> nvent(r) < 0.02), ("open    (>= 0.02)", r -> nvent(r) >= 0.02))
    rs = filter(f, nights); isempty(rs) && (println("  ", lab, ": none"); continue)
    ko = kfit(Cobs, rs); km = kfit(Cmod, rs)
    @printf("  %-32s %3d nights: k obs %.3f, model %.3f (ratio %.2f); R obs %+.1f, model %+.1f ppm/h; mean wind %.2f m/s\n",
            lab, length(rs), ko, km, km / ko, median([Rnight(Cobs, r, ko) for r in rs]) * toppm,
            median([Rnight(Cmod, r, km) for r in rs]) * toppm, mean(nwind.(rs)))
end
rc = filter(r -> nvent(r) < 0.02, nights)
if length(rc) >= 6
    w = nwind.(rc); lo = rc[w .<= median(w)]; hi = rc[w .> median(w)]
    @printf("  closed nights, wind <= median (%.2f m/s): k obs %.3f, model %.3f;   wind > median: k obs %.3f, model %.3f\n",
            median(w), kfit(Cobs, lo), kfit(Cmod, lo), kfit(Cobs, hi), kfit(Cmod, hi))
end
good = filter(r -> nrange(Cobs, r) >= 40.0, nights)    # nights with enough decay to fit k alone (>= ~22 ppm range)
if length(good) >= 5
    ko = [k1(Cobs, r)[1] for r in good]; km = [k1(Cmod, r)[1] for r in good]
    @printf("\n  3. Per-night fits (%d nights with a measured range >= 22 ppm): median k obs %.3f, model %.3f; corr(k model, k obs) %.2f\n",
            length(good), median(ko), median(km), cor(km, ko))
    @printf("     corr(k obs, night vent) %.2f, corr(k obs, wind) %.2f  |  corr(k model, night vent) %.2f, corr(k model, wind) %.2f\n",
            cor(ko, nvent.(good)), cor(ko, nwind.(good)), cor(km, nvent.(good)), cor(km, nwind.(good)))
end
@printf("\n  model vent (night, lamps off) mean %.3f vs logged %.3f\n", mean(Vmod[use]), mean(Vlog[use]))
println("""
  Read: k is the fraction of the excess CO2 exchanged with outside per hour, fitted from how fast each night's CO2
  relaxes -- not from its level, so the daytime dosing error does not enter. k model > k obs = the model exchanges
  too much air at night (leakage or night venting too strong); < = too little. Closed-vent nights isolate leakage.
  R is the net night source (respiration of crop and soil); R model << R obs = respiration too low in the model.""")
