# =============================================================================
# forecast_value.jl -- how much of the SAC->ceiling gap is the value of the
# weather forecast? Evaluate the ALREADY-TRAINED agent (no retraining) across
# forecast skill 0 -> 1 on the TEST years. The agent was trained at skill=1.0,
# so it has forecast inputs in its observation; the campaign TEST used skill=0.
#
#   Run:  julia -t auto --project=. test/forecast_value.jl   (~1-2 min)
# =============================================================================
include(joinpath(@__DIR__, "train_sac.jl"))   # params, gp, bank, actor_sample, DENSITY, Serialization
using Random, Printf, Statistics, Serialization

const TEST_YEARS = [1990, 1997, 2003, 2009, 2014, 2015]
const SKILLS     = [0.0, 0.25, 0.5, 0.75, 1.0]

function rollout(actor, env, year)
    obs = env_reset!(env; density_sched = (_ -> DENSITY), year = year)
    tot = 0.0; done = false
    while !done
        a = vec(Float64.(actor_sample(actor, Float32.(reshape(obs, :, 1)); deterministic = true)[1]))
        obs, r, done, _ = env_step!(env, a); tot += r
    end
    return tot
end

function main()
    afile = isempty(ARGS) ? "sac_actor_best.jls" : ARGS[1]   # e.g. sac_actor_mpc_best.jls
    actor = deserialize(joinpath(@__DIR__, afile))
    @printf("actor: %s\n", afile)
    @printf("forecast value of the trained agent (deployed seed 104), TEST years\n")
    @printf("references: CEM-constant 0.12 | SAC blind (campaign) 1.66 | MPC perfect-foresight ceiling 3.70\n\n")
    @printf("%12s %10s\n", "fcst_skill", "mean TEST")
    for sk in SKILLS
        env = TwinEnv(params, gp, bank; forecast_skill = sk, rng = MersenneTwister(1), crop = :cohort)
        m = mean(rollout(actor, env, y) for y in TEST_YEARS)
        @printf("%12.2f %10.2f\n", sk, m); flush(stdout)
    end
    println("\nRise 0->1 = forecast value the agent can already capture (no retraining).")
    println("Residual below 3.70 at skill=1.0 = gap to the open-loop perfect-foresight optimum.")
end

main()
