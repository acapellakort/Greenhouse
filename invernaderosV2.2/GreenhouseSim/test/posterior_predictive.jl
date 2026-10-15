# =============================================================================
# posterior_predictive.jl -- does the posterior spread cover the data? (Note 64; forward ensemble, no MCMC)
# =============================================================================
# Run from the GreenhouseSim folder (about 1.5-3 s per draw):
#     PRIORS=n54 K=100 julia --project=. test/posterior_predictive.jl
#     PRIORS=n61 K=100 CROP=1 julia --project=. test/posterior_predictive.jl
#
# K parameter draws, evenly spaced over the pooled post-burn chains (burn-in 30 000), each simulated over the full
# season (CL or OL as the PRIORS structure). For every hourly bin:
#   parameter band   5-95 % of the K model hourly means (parameter uncertainty only)
#   predictive band  same plus the fitted observation noise (sigma_T, sigma_RH) drawn per bin
# Coverage = share of measured hourly means inside the band, on TRAIN, VAL, TEST and COMMON-TEST days, day / night.
# A calibrated model gives about 90 % for the predictive band. Coverage well below 90 % means the spread misses
# structural error (bound-pinned parameters, missing physics) or the errors are autocorrelated, not white noise.
# CROP=1: each draw's air T also drives the Note 21 crop model (measured CO2 and light): spread of nodes and yield.
# Compared T includes the measurement terms of the model (n61: s_rad); RH is the air node.
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n61")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
using Random
Random.seed!(2026)
const K    = parse(Int, get(ENV, "K", "100"))
const CROP = get(ENV, "CROP", "0") == "1"

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2"); fc = replace(f, ".jld2" => "_ckpt.jld2")
    src = isfile(f) ? f : isfile(fc) ? fc : ""; isempty(src) && continue
    ch = JLD2.load(src)["chain"]; size(ch, 1) > 31000 || continue
    push!(parts, ch[30001:end, :]); println("seed $s: ", basename(src), "  ", size(ch, 1) - 30000, " post-burn steps")
end
isempty(parts) && error("no $(PRIORS) chains in $POSTDIR")
post = reduce(vcat, parts)
idx  = round.(Int, range(1, size(post, 1), length = K))
kσT  = findfirst(==("sigma_T"), CQ); kσR = findfirst(==("sigma_RH"), CQ)

# hourly bins over the season
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
keys_ = sort(collect(keys(HB)))
okT = .!isnan.(Tobs_g); okT[1] = false; okR = okT .& .!isnan.(RHobs_g)
binmean(v, cols, m) = (k = [c for c in cols if m[c]]; length(k) >= 6 ? mean(v[k]) : NaN)
oT = [binmean(Tobs_g, HB[k], okT) for k in keys_]; oR = [binmean(RHobs_g, HB[k], okR) for k in keys_]
nightbin = [mean(nite[HB[k]]) > 0.5 for k in keys_]; daybin = [mean(dayt[HB[k]]) > 0.5 for k in keys_]
dayof = [k[1] for k in keys_]

MT = fill(NaN, K, length(keys_)); MR = fill(NaN, K, length(keys_)); σT = zeros(K); σR = zeros(K)
nodes = Float64[]; ylds = Float64[]
if CROP
    c21 = JLD2.load(N21_CHAIN); m21 = vec(mean(c21["chain"][c21["burn_in"]:end, 1:end-1], dims = 1))
    gp21 = update_params(gp_base, QoI_G, [m21[findfirst(==(q), c21["QoI"])] for q in QoI_G])
end
t0w = time(); nfail = 0
for (i, r) in enumerate(idx)
    z = post[r, :]
    Tc, RHc, _ = try full_eval(z[1:NM]) catch; (nothing, nothing, nothing) end
    if Tc === nothing; global nfail += 1; continue; end
    MT[i, :] = [binmean(Tc, HB[k], okT) for k in keys_]; MR[i, :] = [binmean(RHc, HB[k], okR) for k in keys_]
    σT[i] = z[kσT]; σR[i] = z[kσR]
    if CROP     # crop with this draw's air T (s_rad excluded: the crop sees the air, not the sensor)
        us_ = simdays(update_params(params_cl, CQ_MODEL, z[1:NM]), days, u0_cl, RELTOL); Tair = us_[2, :]
        for d in 0:(N_DAYS - 1); M = work[d + 1]
            for k in axes(M, 2); v = Tair[samp_col[d + 1][k]]; M[1, k] = isnan(v) ? samp_meas[d + 1][1, k] : v; end
        end
        rr = simulate_growth_measured(params, gp21, season_B; LAI0 = 0.5); push!(nodes, rr.node[end]); push!(ylds, rr.yield_FW[end])
    end
    (i % 10 == 0) && @printf("  draw %d/%d  %.0f s\n", i, K, time() - t0w)
end
ok = .!isnan.(MT[:, 1]); @printf("\n%d draws simulated, %d failed\n", count(ok), nfail)
MT = MT[ok, :]; MR = MR[ok, :]; σT = σT[ok]; σR = σR[ok]
q(v, p) = (w = filter(!isnan, v); isempty(w) ? NaN : quantile(w, p))
# predictive draws: each ensemble member's bin mean plus N(0, sigma) of that member
PT = MT .+ σT .* randn(size(MT)); PR = MR .+ σR .* randn(size(MR))
loT = [q(MT[:, j], 0.05) for j in axes(MT, 2)]; hiT = [q(MT[:, j], 0.95) for j in axes(MT, 2)]
plT = [q(PT[:, j], 0.05) for j in axes(PT, 2)]; phT = [q(PT[:, j], 0.95) for j in axes(PT, 2)]
plR = [q(PR[:, j], 0.05) for j in axes(PR, 2)]; phR = [q(PR[:, j], 0.95) for j in axes(PR, 2)]
loR = [q(MR[:, j], 0.05) for j in axes(MR, 2)]; hiR = [q(MR[:, j], 0.95) for j in axes(MR, 2)]
cov_(o, lo, hi, m) = (k = m .& .!isnan.(o) .& .!isnan.(lo); count(k) == 0 ? NaN : 100 * mean((o[k] .>= lo[k]) .& (o[k] .<= hi[k])))
wid(lo, hi, m) = (k = m .& .!isnan.(lo); mean(hi[k] .- lo[k]))
sets = [("TRAIN", TRAIN_DAYS), ("VAL", VAL_DAYS), ("TEST", TEST_DAYS), ("COMMON TEST", COMMON_TEST)]
println("\n", "="^120)
@printf("  %s posterior predictive, %d draws.  Coverage (%%) of measured hourly means by the 90%% band; mean band width\n", PRIORS, count(ok))
@printf("  %-12s | %-34s | %-34s\n", "", "air T: parameter band / predictive band", "RH: parameter band / predictive band")
for (nm_, ds) in sets
    m = in.(dayof, Ref(Set(ds)))
    @printf("  %-12s | %5.1f %% (%.2f K) / %5.1f %% (%.2f K) | %5.1f %% (%.1f %%) / %5.1f %% (%.1f %%)\n", nm_,
            cov_(oT, loT, hiT, m), wid(loT, hiT, m), cov_(oT, plT, phT, m), wid(plT, phT, m),
            cov_(oR, loR, hiR, m), wid(loR, hiR, m), cov_(oR, plR, phR, m), wid(plR, phR, m))
    @printf("  %-12s | day %5.1f %%  night %5.1f %% (predictive T)          | day %5.1f %%  night %5.1f %% (predictive RH)\n", "",
            cov_(oT, plT, phT, m .& daybin), cov_(oT, plT, phT, m .& nightbin), cov_(oR, plR, phR, m .& daybin), cov_(oR, plR, phR, m .& nightbin))
end
# autocorrelation of the residual of the posterior-mean bins (white noise -> ~0)
mT = vec(mean(MT, dims = 1)); res = oT .- mT; kk = .!isnan.(res)
r1 = cor(res[kk][1:end-1], res[kk][2:end]); r24 = cor(res[kk][1:end-24], res[kk][25:end])
@printf("\n  residual of the ensemble-mean air T: lag-1 h autocorrelation %.2f, lag-24 h %.2f (white noise: ~0)\n", r1, r24)
if CROP && !isempty(nodes)
    @printf("  crop spread over %d draws: nodes %.1f [%.1f, %.1f], final yield %.2f [%.2f, %.2f] kg m-2 (5-95 %%)\n", length(nodes),
            median(nodes), quantile(nodes, 0.05), quantile(nodes, 0.95), median(ylds), quantile(ylds, 0.05), quantile(ylds, 0.95))
end
println("""
  Read: predictive coverage near 90 % = the posterior plus noise is an honest uncertainty for uses B-E.
  Much lower = add structural uncertainty (spread across model versions, or an autocorrelated error term).
  The parameter band alone is expected to be narrow; its width is what parameter uncertainty contributes.""")
