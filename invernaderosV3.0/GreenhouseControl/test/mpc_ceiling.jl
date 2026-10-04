# =============================================================================
# mpc_ceiling.jl -- perfect-foresight open-loop ceiling for the twin.
#
# The twin is deterministic given a weather year, so the best OPEN-LOOP action
# schedule with the weather REVEALED is the per-year optimum -- an upper bound on
# any causal controller. We approximate it per TEST year with CEM over a block
# schedule of the 11 setpoints, WARM-STARTED from the trained SAC agent's own
# realized trajectory (so the search starts at a strong point and can only add
# the value of perfect foresight + reoptimization).
#
# Reports per year: CEM-constant (WARM), SAC closed-loop, SAC-block-avg (a
# resolution check), and the MPC ceiling. Gap SAC->ceiling = headroom left.
#
#   Run:  julia -t auto --project=. test/mpc_ceiling.jl
# =============================================================================
include(joinpath(@__DIR__, "train_sac.jl"))     # params, gp, bank, actor_sample, ACT_*, DENSITY, WARM_A, Serialization
using Random, Printf, Statistics, Serialization

const SMOKE = "smoke" in ARGS   # fast wiring check: 2 years, tiny CEM (~2 min)
const TEST_YEARS = SMOKE ? [1990, 2014] : [1990, 1997, 2003, 2009, 2014, 2015]
const NB    = 10                       # schedule blocks over the season (10-day blocks at ndays=100)
const POP   = SMOKE ? 8  : 48
const NEL   = SMOKE ? 3  : 8
const ITERS = SMOKE ? 2  : 20
const STD0  = 0.10

day_block(day, ndays, nb) = min(nb, 1 + (day * nb) ÷ ndays)   # day 0-based -> block 1..nb

"Apply a per-day action matrix A (ACT_DIM x ndays) to a fixed year; return total profit."
function eval_daily(env, A, year)
    env_reset!(env; density_sched = (_ -> DENSITY), year = year)
    tot = 0.0; done = false; day = 0
    while !done
        _, r, done, _ = env_step!(env, clamp.(A[:, day+1], 0.0, 1.0)); tot += r; day += 1
    end
    return tot
end

expand(S, ndays) = hcat((S[:, day_block(d, ndays, NB)] for d in 0:ndays-1)...)  # (ACT_DIM x ndays)
eval_block(env, z, year, ndays) = eval_daily(env, expand(reshape(z, ACT_DIM, NB), ndays), year)

"SAC closed-loop rollout: (profit, realized daily actions ACT_DIM x ndays)."
function sac_rollout(env, actor, year)
    obs = env_reset!(env; density_sched = (_ -> DENSITY), year = year)
    tot = 0.0; done = false; acts = Vector{Vector{Float64}}()
    while !done
        a = vec(Float64.(actor_sample(actor, Float32.(reshape(obs, :, 1)); deterministic = true)[1]))
        push!(acts, a); obs, r, done, _ = env_step!(env, a); tot += r
    end
    return tot, reduce(hcat, acts)
end

"Block-average a daily action matrix into a flat NB*ACT_DIM schedule."
function block_init(A, ndays)
    M = Matrix{Float64}(undef, ACT_DIM, NB)
    for b in 1:NB
        cols = [d+1 for d in 0:ndays-1 if day_block(d, ndays, NB) == b]
        M[:, b] = isempty(cols) ? A[:, end] : vec(mean(A[:, cols], dims = 2))
    end
    return vec(M)
end

"Warm-started diagonal CEM over the block schedule for one fixed (revealed) year."
function cem_year(env, year, z0, ndays)
    dim = length(z0); μ = copy(z0); σ = fill(STD0, dim)
    best_z = copy(z0); best_v = eval_block(env, z0, year, ndays)
    for _ in 1:ITERS
        Z = [clamp.(μ .+ σ .* randn(dim), 0.0, 1.0) for _ in 1:POP]
        V = [eval_block(env, z, year, ndays) for z in Z]
        idx = partialsortperm(V, 1:NEL, rev = true)
        E = reduce(hcat, Z[idx])
        μ = vec(mean(E, dims = 2)); σ = vec(std(E, dims = 2)) .+ 1e-3
        if V[idx[1]] > best_v; best_v = V[idx[1]]; best_z = Z[idx[1]]; end
    end
    return best_v
end

function main()
    actor = deserialize(joinpath(@__DIR__, "sac_actor_best.jls"))
    env   = TwinEnv(params, gp, bank; forecast_skill = 0.0, rng = MersenneTwister(1), crop = :cohort)
    ndays = env.ndays
    zc = vec(repeat(reshape(Float64.(WARM_A), ACT_DIM, 1), 1, NB))   # constant CEM program, all blocks

    @printf("%6s %10s %10s %12s %10s %10s\n", "year", "CEM-const", "SAC", "SAC-blockavg", "MPC-ceil", "gap SAC->ceil")
    rows = NTuple{5,Float64}[]
    for y in TEST_YEARS
        cem_c = eval_block(env, zc, y, ndays)
        sac_v, A = sac_rollout(env, actor, y)
        z0 = block_init(A, ndays)
        sba = eval_block(env, z0, y, ndays)          # SAC block-averaged (resolution check)
        ceil = cem_year(env, y, z0, ndays)
        push!(rows, (cem_c, sac_v, sba, ceil, ceil - sac_v))
        @printf("%6d %10.2f %10.2f %12.2f %10.2f %10.2f\n", y, cem_c, sac_v, sba, ceil, ceil - sac_v); flush(stdout)
    end
    m = [mean(getindex.(rows, i)) for i in 1:5]
    @printf("%6s %10.2f %10.2f %12.2f %10.2f %10.2f\n", "mean", m...)
    println("\nNotes: SAC-blockavg near SAC => 10-day blocks retain the agent's detail (ceiling is fair).")
    println("gap SAC->ceil small => SAC near the perfect-foresight frontier; large => headroom remains.")
end

main()
