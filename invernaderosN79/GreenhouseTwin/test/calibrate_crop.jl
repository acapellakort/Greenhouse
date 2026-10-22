# =============================================================================
# calibrate_crop.jl -- crop calibration (n105 set-up), flag-free version of V2.2 calibrate_canopy_growth.jl PARBELOW=1
# =============================================================================
# Run from invernaderosN79/GreenhouseTwin:
#   MODE=centre julia --project=. test/calibrate_crop.jl                                   (seconds: fit at the prior centre)
#   MODE=smoke  julia --project=. test/calibrate_crop.jl                                   (2000 steps)
#   MODE=full SEED=1 CHAIN=500000 BURN=150000 julia --project=. test/calibrate_crop.jl     (long: ask first)
# Driver: measured AiCU air T and CO2, canopy-top PAR by the organisers' formula (transmission 0.65, lamps 187 umol),
# times the both-faces factor g(LAI). Observations: weekly yield, weekly harvested fruit number, leaves per stem,
# fruit development time. Free parameters and priors: configfiles/crop_n105.json.
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin
using JTwalk, JLD2, Printf, Random, Statistics
const GT = GreenhouseTwin

const MODE  = get(ENV, "MODE", "centre"); MODE in ("centre", "smoke", "full") || error("MODE = centre | smoke | full")
const SEED  = parse(Int, get(ENV, "SEED", "1")); Random.seed!(SEED)
const CHAIN = parse(Int, get(ENV, "CHAIN", MODE == "smoke" ? "2000" : "500000"))
const BURN  = parse(Int, get(ENV, "BURN", "150000"))
const INV   = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CJ    = joinpath(@__DIR__, "..", "configfiles", "climate_n79.json")
const KJ    = joinpath(@__DIR__, "..", "configfiles", "crop_n105.json")
const pc  = load_params(CJ)
const gp0 = load_params(KJ)
const PR  = priors_of(KJ)
const Q   = vcat(GT.CROP_QOI, GT.CROP_NOISE)
sort(inferred_names(KJ)) == sort(Q) || error("crop_n105.json inferred set differs from CROP_QOI + noise")
const O   = crop_observations(INV; transmission = gp0.transmission, lamp_umol = gp0.lamp_umol)
const centre = [PR[q][2] for q in Q]
ix(q) = findfirst(==(q), Q)
supp(z) = all(PR[Q[k]][1] < z[k] < PR[Q[k]][3] for k in eachindex(z)) && z[ix("set_r_high")] > z[ix("set_r_low")] + 0.05
energy(z) = crop_energy(z, pc, gp0, O)

function report(z, tag)
    r = GT.crop_run(z, pc, gp0, O); ry, rl, rf, rn = GT.crop_residuals(r, O)
    nh = sum(r.harv_n)
    @printf("\n  %s: final yield %.2f kg FW m-2 (obs %.2f) | weekly-yield RMSE %.3f | leaves RMSE %.2f | fruit-age RMSE %.2f d\n",
            tag, r.yield_FW[end], GT._cum_at(O.days_obs, O.y_obs, O.last), sqrt(mean(abs2, ry)), sqrt(mean(abs2, rl)), sqrt(mean(abs2, rf)))
    @printf("  fruits harvested %.0f m-2 (obs A+B %.0f), mean fruit %.0f g FW; LAI max %.2f end %.2f; retained leaves end %.0f\n",
            nh, sum(O.dn), nh > 0 ? 1000 * r.yield_FW[end] / nh : NaN, maximum(r.LAI), r.LAI[end], r.leaves[end])
end

println("\nCrop calibration (GreenhouseTwin, n105 set-up) -- mode $MODE seed $SEED")
report(centre, "prior centre"); @printf("  energy at the centre %.2f\n", energy(centre))
MODE == "centre" && exit()
jit() = (for _ in 1:1000; x = centre .* (1 .+ 0.05 .* randn(length(Q))); supp(x) && return x; end; centre)
obj = jtwalk(n = length(Q), U = energy, Supp = supp)
t0 = time(); Run!(obj, T = CHAIN, x0 = jit(), xp0 = jit())
@printf("\n%d steps in %.0f s\n", CHAIN, time() - t0)
POST = joinpath(@__DIR__, "posteriors"); mkpath(POST)
fn = joinpath(POST, "AiCU_crop_n105_twin_seed$(SEED)$(MODE == "smoke" ? "_smoke" : "").jld2")
JLD2.@save fn chain = obj.Output QoI = Q burn_in = BURN
println("saved ", fn)
b = MODE == "smoke" ? max(1, CHAIN ÷ 2) : BURN
post = obj.Output[b:end, 1:end-1]; pm = vec(mean(post, dims = 1)); ps = vec(std(post, dims = 1))
println("\nposterior (mean ± sd):"); for (k, q) in enumerate(Q); @printf("  %-16s %10.4g ± %.3g\n", q, pm[k], ps[k]); end
report(pm, "posterior mean")
