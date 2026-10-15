# =============================================================================
# residual_error_model.jl -- an honest error model for forecasts: day/night sigma + autocorrelation (Note 65)
# =============================================================================
# Run from the GreenhouseSim folder (~1-2 min per model, forward only, no MCMC):
#     PRIORS=n54 julia --project=. test/residual_error_model.jl
#     PRIORS=n61 julia --project=. test/residual_error_model.jl
#
# Note 64: the posterior adds almost no width (1-2 % of the variance); the residual has lag-1 h autocorrelation
# ~0.8 and lag-24 h ~0.47; one sigma_T over-covers nights and under-covers days. Here, from the posterior-mean
# run (hourly means, residual e = measured - model, compared T as in the likelihood):
#   fit on TRAIN only:   sigma_day, sigma_night (RMS incl. bias), and the autocorrelation of e / sigma at lags
#                        1..72 h, matched by  M2 AR(1)  and  M3 fast + slow AR (acf = w rf^k + (1-w) rs^k)
#   test on all sets:    coverage of 90 % bands and the error-to-predicted-SD ratio (target 1.0) for the hourly
#                        value, the 24-h mean and the 72-h mean, for
#                        M0 posterior sigma, independent  |  M1 day/night sigma, independent
#                        M2 day/night sigma, AR(1)        |  M3 day/night sigma, fast + slow AR
# Same for RH. The printed M3 line is the error model to add to a forecast (use B) or to RL/MPC rollouts (use C).
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n54")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
using Random
Random.seed!(64)
const NSIM = parse(Int, get(ENV, "NSIM", "2000"))

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2"); fc = replace(f, ".jld2" => "_ckpt.jld2")
    src = isfile(f) ? f : isfile(fc) ? fc : ""; isempty(src) && continue
    ch = JLD2.load(src)["chain"]; size(ch, 1) > 31000 || continue
    push!(parts, ch[30001:end, 1:length(CQ)]); println("seed $s: ", basename(src))
end
isempty(parts) && error("no $(PRIORS) chains in $POSTDIR")
pmean = vec(mean(reduce(vcat, parts), dims = 1))
sigT_post = pmean[findfirst(==("sigma_T"), CQ)]; sigR_post = pmean[findfirst(==("sigma_RH"), CQ)]
Tcmp, RHm, _ = full_eval(pmean[1:NM])

# ── hourly series on an absolute hour index (h = 0 .. NH-1; day = h ÷ 24) ─────────
hix = [floor(Int, (predtimes[c] - t0) / 3600) for c in 1:NCOL]
NH  = maximum(hix[2:end]) + 1
colsof = [Int[] for _ in 1:NH]
for c in 2:NCOL; push!(colsof[hix[c] + 1], c); end
okTc = .!isnan.(Tobs_g); okTc[1] = false; okRc = okTc .& .!isnan.(RHobs_g)
function hseries(mv, ov, ok)
    m = fill(NaN, NH); o = fill(NaN, NH)
    for h in 1:NH
        cs = colsof[h]; isempty(cs) && continue; m[h] = mean(mv[cs])
        k = [c for c in cs if ok[c]]; length(k) >= 6 && (o[h] = mean(ov[k]))
    end
    m, o
end
isday_h = [isempty(colsof[h]) ? false : mean(dayt[colsof[h]]) > 0.5 for h in 1:NH]
day_h   = [(h - 1) ÷ 24 for h in 1:NH]
inds(ds) = in.(day_h, Ref(Set(ds)))
SETS = [("TRAIN", TRAIN_DAYS), ("VAL", VAL_DAYS), ("TEST", TEST_DAYS), ("COMMON", COMMON_TEST)]

function acf_u(z, mask, k)       # uncentred autocorrelation (the bias is part of the error), pairs inside mask
    a = 0.0; b = 0.0; c = 0.0; n = 0
    for t in 1:(length(z) - k)
        (mask[t] && mask[t + k] && !isnan(z[t]) && !isnan(z[t + k])) || continue
        a += z[t] * z[t + k]; b += z[t]^2; c += z[t + k]^2; n += 1
    end
    n < 20 ? NaN : a / sqrt(b * c)
end

function arpaths(r, n)           # NSIM x n stationary AR(1) paths with unit variance
    X = Matrix{Float64}(undef, NSIM, n); X[:, 1] = randn(NSIM); s = sqrt(1 - r^2)
    for t in 2:n
        X[:, t] = r .* X[:, t - 1] .+ s .* randn(NSIM)
    end
    return X
end
function simz(kind, p, n)        # NSIM x n unit-variance stationary paths
    kind == :iid && return randn(NSIM, n)
    kind == :ar1 && return arpaths(p.r, n)
    return sqrt(p.w) .* arpaths(p.rf, n) .+ sqrt(1 - p.w) .* arpaths(p.rs, n)
end

function windows(ds, w)          # non-overlapping w-hour windows, starting at midnight, inside the day set
    out = Vector{UnitRange{Int}}(); s = Set(ds)
    for d in sort(collect(ds))
        (w > 24 && (d - minimum(ds)) % (w ÷ 24) != 0) && continue
        r = (24d + 1):(24d + w); last(r) <= NH || continue
        all(in(s), unique(day_h[r])) && push!(out, r)
    end
    w == 1 ? [h:h for d in ds for h in (24d + 1):min(24d + 24, NH)] : out
end

function evaluate(name, m, o, sigpost, unit)
    e = o .- m; tr = inds(TRAIN_DAYS)
    sd = sqrt(mean(filter(!isnan, e[tr .& isday_h]) .^ 2)); sn = sqrt(mean(filter(!isnan, e[tr .& .!isday_h]) .^ 2))
    sig = [isday_h[h] ? sd : sn for h in 1:NH]
    z = e ./ sig
    lags = 1:72; ac = [acf_u(z, tr, k) for k in lags]; okl = .!isnan.(ac)
    r1 = ac[1]
    best = (Inf, 0.0, 0.0, 0.0)
    for rf in 0.30:0.01:0.98, rs in 0.950:0.001:0.999, w in 0.0:0.02:1.0
        f = sum(abs2, (w .* rf .^ lags .+ (1 - w) .* rs .^ lags .- ac)[okl]); f < best[1] && (best = (f, rf, rs, w))
    end
    _, rf, rs, w = best
    println("\n", "-"^112)
    @printf("  %s  (%s)   TRAIN fit:  sigma day %.2f, night %.2f %s  (posterior sigma %.2f);  TRAIN bias day %+.2f, night %+.2f\n",
            name, PRIORS, sd, sn, unit, sigpost, mean(filter(!isnan, e[tr .& isday_h])), mean(filter(!isnan, e[tr .& .!isday_h])))
    @printf("  M2 AR(1): rho %.3f      M3: w %.2f x fast rho %.3f (%.1f h) + %.2f x slow rho %.4f (%.1f days)\n",
            r1, w, rf, -1 / log(rf), 1 - w, rs, -1 / log(rs) / 24)
    @printf("  autocorrelation of e/sigma  %-8s", "lag (h):"); for k in (1, 2, 3, 6, 12, 24, 48, 72); @printf("%7d", k); end; println()
    for (lab, msk) in (("TRAIN", tr), ("VAL", inds(VAL_DAYS)), ("TEST", inds(TEST_DAYS)))
        @printf("  %-36s", "  measured " * lab); for k in (1, 2, 3, 6, 12, 24, 48, 72); @printf("%7.2f", acf_u(z, msk, k)); end; println()
    end
    @printf("  %-36s", "  M2 AR(1)"); for k in (1, 2, 3, 6, 12, 24, 48, 72); @printf("%7.2f", r1^k); end; println()
    @printf("  %-36s", "  M3 fast + slow"); for k in (1, 2, 3, 6, 12, 24, 48, 72); @printf("%7.2f", w * rf^k + (1 - w) * rs^k); end; println()

    MODELS = [("M0 posterior sigma, indep.", fill(sigpost, NH), :iid, nothing),
              ("M1 day/night sigma, indep.", sig, :iid, nothing),
              ("M2 day/night sigma, AR(1)", sig, :ar1, (r = r1,)),
              ("M3 day/night, fast+slow AR", sig, :two, (rf = rf, rs = rs, w = w))]
    @printf("\n  coverage of the 90 %% band (%%) / ratio RMS(error) : predicted SD (target 1.0); n = windows\n")
    @printf("  %-28s %-7s", "", "window"); for (s, _) in SETS; @printf(" | %-20s", s); end; println()
    for (lab, sg, kind, p) in MODELS
        E = simz(kind, p, NH) .* sg'
        for w_ in (1, 24, 72)
            @printf("  %-28s %-7s", w_ == 1 ? lab : "", w_ == 1 ? "1 h" : "$(w_) h")
            for (_, ds) in SETS
                cov = 0; rat = 0.0; n = 0
                for r in windows(ds, w_)
                    k = [h for h in r if !isnan(o[h])]; length(k) >= 0.8 * length(r) || continue
                    eb = mean(e[k]); sim = vec(mean(E[:, k], dims = 2))
                    lo, hi = quantile(sim, 0.05), quantile(sim, 0.95)
                    cov += (lo <= eb <= hi); rat += (eb / std(sim))^2; n += 1
                end
                if n == 0; @printf(" | %-20s", "-"); else; @printf(" | %5.1f %% / %4.2f n%-4d", 100cov / n, sqrt(rat / n), n); end
            end
            println()
        end
    end
    @printf("\n  => error model for forecasts / rollouts (%s, %s):  e = s(t) * (sqrt(%.2f) a + sqrt(%.2f) b),  s = %.2f day / %.2f night %s,\n",
            name, PRIORS, w, 1 - w, sd, sn, unit)
    @printf("     a: AR(1) rho %.3f per hour,  b: AR(1) rho %.4f per hour  (stationary, unit variance)\n", rf, rs)
end

mT, oT_ = hseries(Tcmp, Tobs_g, okTc); mR, oR_ = hseries(RHm, RHobs_g, okRc)
println("\n", "="^112)
@printf("  %s residual error model -- posterior mean, %s, %d simulated error paths\n", PRIORS,
        OLP ? "OPEN loop (logged vents)" : "CLOSED loop (model vents)", NSIM)
evaluate("air T", mT, oT_, sigT_post, "K")
evaluate("RH", mR, oR_, sigR_post, "%")
println("""
  Read: M0 is what Note 64 used. A model is adequate for use B/C when coverage is ~90 % and the ratio ~1.0 on
  VAL / TEST at ALL three windows -- the 24-h and 72-h rows are what crop development, energy and RL returns see.
  Ratio > 1 = the model is overconfident at that window (errors persist more than it assumes).""")
