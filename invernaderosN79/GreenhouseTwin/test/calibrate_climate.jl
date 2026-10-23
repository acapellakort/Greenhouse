# =============================================================================
# calibrate_climate.jl -- n79 climate calibration (t-walk), flag-free version of V2.2 run_climate_calibration_cl.jl
# =============================================================================
# Run from invernaderosN79/GreenhouseTwin:
#   CROPTAG=n105 RELTOL=1e-4 julia --project=. test/calibrate_climate.jl smoke 2>&1 | tee climcal_smoke.log   # ALWAYS FIRST
#   CROPTAG=n105 SEED=1 CHAIN=100000 BURN=40000 RELTOL=1e-4 \
#       julia --project=. test/calibrate_climate.jl full 2>&1 | tee climcal_seed1.log                       # long: ask first
# Free parameters, priors and fixed values: configfiles/climate_n79.json (status / prior). Likelihood: hourly air T
# and RH, dark-hour CO2 on the TRAIN blocks (src/climate_calibration.jl). Crop inputs: test/posteriors/AiCU_crop_series_<CROPTAG>.csv
# (written by test/crop_series.jl; falls back to the V2.2 posteriors folder).
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin
using JTwalk, JLD2, Printf, Random, Statistics, DataFrames
const GT = GreenhouseTwin

const MODE    = isempty(ARGS) ? "smoke" : ARGS[1]
MODE in ("smoke", "full") || error("usage: calibrate_climate.jl [smoke|full]")
const SEED    = parse(Int, get(ENV, "SEED", "1")); Random.seed!(SEED)
const CHAIN   = parse(Int, get(ENV, "CHAIN", "100000"))
const BURN    = parse(Int, get(ENV, "BURN", "40000"))
const RELTOL  = parse(Float64, get(ENV, "RELTOL", "1e-4"))
const CROPTAG  = get(ENV, "CROPTAG", "n105")
const CLIMTAG  = get(ENV, "CLIMTAG", "n79")
const NSMOKE  = 200
const INV     = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CJ      = joinpath(@__DIR__, "..", "configfiles", "climate_n79.json")
const POST    = joinpath(@__DIR__, "posteriors"); mkpath(POST)
cs_new = joinpath(POST, "AiCU_crop_series_$(CROPTAG).csv")
cs_old = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "test", "posteriors", "AiCU_crop_series_$(CROPTAG).csv")
const CROPCSV = isfile(cs_new) ? cs_new : cs_old

const p0  = load_params(CJ)
const PR  = priors_of(CJ)
const Q   = vcat(GT.CLIMATE_QOI, GT.CLIMATE_NOISE)
sort(inferred_names(CJ)) == sort(Q) || error("climate_n79.json inferred set differs from CLIMATE_QOI + noise")
const S   = load_agc_season(INV; crop_series = CROPCSV)
const BL  = climate_blocks(S)
const centre = [PR[q][2] for q in Q]
supp(z) = all(PR[Q[k]][1] < z[k] < PR[Q[k]][3] for k in eachindex(z))
energy(z; rtol = RELTOL) = climate_energy(z, p0, S, BL; rtol = rtol)
function jitter(c)
    for _ in 1:1000
        x = c .* (1 .+ 0.05 .* randn(length(c))); supp(x) && return x
    end
    error("no start point in the support")
end

println("="^78)
@printf("n79 climate calibration (GreenhouseTwin) -- mode %s, seed %d, reltol %.0e, crop series %s\n", MODE, SEED, RELTOL, CROPCSV)
@printf("TRAIN bins: T %d, RH %d, dark CO2 %d, day VP %d, day CO2 %d (%d simulated days per energy)\n",
        sum(length(b.obsT) for b in BL), sum(length(b.obsR) for b in BL), sum(length(b.obsC) for b in BL),
        sum(length(b.binsVP) for b in BL), sum(length(b.binsCd) for b in BL), sum(length(b.days) for b in BL))
println("="^78)

if MODE == "smoke"
    xa = copy(centre); Random.seed!(SEED + 100); xb = jitter(centre); Random.seed!(SEED)
    energy(xa; rtol = 1e-3)
    println("\n[smoke 1] tolerance: energy at the prior centre xa and a 5%-jittered point xb")
    @printf("  %-8s %12s %12s %12s %10s\n", "reltol", "U(xa)", "U(xb)", "U(xb)-U(xa)", "ms/eval")
    res = Dict{Float64,Tuple{Float64,Float64}}()
    for rt in (1e-3, 1e-4, 1e-5)
        te = @elapsed (Ua = energy(xa; rtol = rt)); Ub = energy(xb; rtol = rt); res[rt] = (Ua, Ub)
        @printf("  %-8.0e %12.2f %12.2f %12.2f %10.0f\n", rt, Ua, Ub, Ub - Ua, 1e3te)
    end
    d3 = (res[1e-3][2] - res[1e-3][1]) - (res[1e-5][2] - res[1e-5][1]); d4 = (res[1e-4][2] - res[1e-4][1]) - (res[1e-5][2] - res[1e-5][1])
    @printf("  error of the energy difference vs 1e-5: 1e-3 %+.2f, 1e-4 %+.2f -> %s\n", d3, d4,
            abs(d3) < 1 ? "1e-3 adequate" : abs(d4) < 1 ? "use RELTOL=1e-4" : "use RELTOL=1e-5")
    obj = jtwalk(n = length(Q), U = energy, Supp = supp)
    tw = @elapsed Run!(obj, T = NSMOKE, x0 = jitter(centre), xp0 = jitter(centre))
    O = obj.Output; acc = mean(any(O[k+1, 1:end-1] .!= O[k, 1:end-1]) for k in 1:size(O, 1)-1)
    @printf("\n[smoke 2] %d steps in %.1f s (%.0f ms/step), moved fraction %.2f; projected %d steps: %.1f h\n",
            NSMOKE, tw, 1e3tw / NSMOKE, acc, CHAIN, tw / NSMOKE * CHAIN / 3600)
    _, tab = season_eval(update_params(p0, GT.CLIMATE_QOI, centre[1:length(GT.CLIMATE_QOI)]), S; rtol = RELTOL)
    println("\n[smoke 3] full season at the PRIOR CENTRE (not calibrated); K, hPa, %, mg m-3; bias = model - observed")
    show(stdout, tab; allrows = true, eltypes = false); println()
    println("\nSmoke done. Paste this output before launching `full`.")
    exit(0)
end

obj = jtwalk(n = length(Q), U = energy, Supp = supp)
@time Run!(obj, T = CHAIN, x0 = jitter(centre), xp0 = jitter(centre))
fn = joinpath(POST, "AiCU_climate_$(CLIMTAG)_twin_seed$(SEED).jld2")
JLD2.@save fn chain = obj.Output QoI = Q burn_in = BURN seed = SEED reltol = RELTOL croptag = CROPTAG climtag = CLIMTAG
println("saved ", fn)
post = obj.Output[BURN+1:end, 1:end-1]; pm = vec(mean(post, dims = 1))
println("\nposterior (mean, 95 %):")
for (i, q) in enumerate(Q)
    lo, _, hi = PR[q]; qs = quantile(post[:, i], [0.025, 0.975]); span = hi - lo
    flag = qs[1] - lo < 0.02span ? "AT LOWER BOUND" : hi - qs[2] < 0.02span ? "AT UPPER BOUND" : ""
    @printf("  %-14s %10.4g  [%10.4g, %10.4g]  %s\n", q, pm[i], qs[1], qs[2], flag)
end
_, tab = season_eval(update_params(p0, GT.CLIMATE_QOI, pm[1:length(GT.CLIMATE_QOI)]), S; rtol = RELTOL)
println("\nfull season at the POSTERIOR MEAN:"); show(stdout, tab; allrows = true, eltypes = false); println()
