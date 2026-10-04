# =============================================================================
# ensemble_eval.jl -- action-average the top-k (by VAL) campaign actors and
#   score the ensemble on the held-out TEST years vs CEM and the best single.
#   Run AFTER train_sac_long.jl:   julia --project=. test/ensemble_eval.jl [k]
# =============================================================================
include(joinpath(@__DIR__, "train_sac.jl"))       # machinery + params/gp/bank/DENSITY/ACT_DIM/WARM_A
using Serialization, Statistics, Printf, DelimitedFiles

const TEST = [1990,1997,2003,2009,2014,2015]
K = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5

man   = readdlm(joinpath(@__DIR__, "campaign_manifest.csv"), ',', skipstart=1)
seeds = Int.(man[:,1]); vals = Float64.(man[:,2])
order = sortperm(vals, rev=true)
topk  = order[1:min(K, length(order))]
@printf("ensembling top %d of %d seeds by VAL (seeds %s)\n",
        length(topk), length(seeds), string(seeds[topk]))
actors = [deserialize(joinpath(@__DIR__, "campaign_actors", "seed_$(seeds[i]).jls")) for i in topk]

evenv = TwinEnv(params, gp, bank; crop=:cohort, forecast_skill=0.0)
det_of(a) = obs -> Float64.((tanh.(a(Float32.(reshape(obs,:,1)))[1:ACT_DIM,1]) .+ 1) ./ 2)
ens(obs)  = (as = [det_of(a)(obs) for a in actors]; reduce(+, as) ./ length(as))
function profits(actfn, years)
    map(years) do y
        obs = env_reset!(evenv; density_sched=(_->DENSITY), year=y); tot=0.0; done=false
        while !done; obs, r, done, _ = env_step!(evenv, actfn(obs)); tot += r; end
        tot
    end
end

cem  = profits(obs->Float64.(WARM_A), TEST)
best = profits(det_of(actors[1]), TEST)   # single best-VAL actor
ense = profits(ens, TEST)
@printf("\n%-6s %8s %8s %10s\n", "year", "CEM", "best-1", "ens($(length(actors)))")
for (i,y) in enumerate(TEST); @printf("%-6d %8.2f %8.2f %10.2f\n", y, cem[i], best[i], ense[i]); end
println("-"^36)
@printf("%-6s %8.2f %8.2f %10.2f   (mean TEST EUR/m2)\n", "mean", mean(cem), mean(best), mean(ense))
@printf("%-6s %8s %8.2f %10.2f   (worst TEST year)\n", "worst", "", minimum(best), minimum(ense))
