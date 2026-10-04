# =============================================================================
# gen_mpc_demos.jl -- forecast-driven EXPERT demonstrations for MPC->RL distill.
#
# For each TRAIN year: warm-start from the deployed SAC trace, refine an open-loop
# schedule with the weather revealed (perfect-foresight CEM), then roll it out and
# record (obs, expert_action) pairs. The obs carries a REALISTIC (skill=1) forecast
# while the action used full foresight -- so imitating these teaches the agent to
# map its noisy near-term forecast to anticipatory actions. Saved -> mpc_demos.jls
#
#   smoke:  julia -t auto --project=. test/gen_mpc_demos.jl smoke   (~1 min, 2 yrs)
#   full :  julia -t auto --project=. test/gen_mpc_demos.jl         (~1.5 h, TRAIN)
# =============================================================================
include(joinpath(@__DIR__, "train_sac.jl"))   # params, gp, bank, actor_sample, ACT_*, DENSITY, Serialization
using Random, Printf, Statistics, Serialization

const SMOKE = "smoke" in ARGS
const TEST_YEARS = [1990, 1997, 2003, 2009, 2014, 2015]
const VAL_YEARS  = [1992, 1999, 2005, 2011, 2016, 2017]
TRAIN_YEARS = setdiff(sort(bank.years), vcat(TEST_YEARS, VAL_YEARS))
SMOKE && (TRAIN_YEARS = TRAIN_YEARS[1:2])

const NB    = 10
const POP   = SMOKE ? 8 : 32
const NEL   = SMOKE ? 3 : 6
const ITERS = SMOKE ? 2 : 10
const STD0  = 0.10

day_block(day, ndays, nb) = min(nb, 1 + (day * nb) ÷ ndays)
expand(S, ndays) = hcat((S[:, day_block(d, ndays, NB)] for d in 0:ndays-1)...)

function eval_block(env, z, year, ndays)
    A = expand(reshape(z, ACT_DIM, NB), ndays)
    env_reset!(env; density_sched = (_ -> DENSITY), year = year)
    tot = 0.0; done = false; day = 0
    while !done
        _, r, done, _ = env_step!(env, clamp.(A[:, day+1], 0.0, 1.0)); tot += r; day += 1
    end
    return tot
end

function sac_trace(env, actor, year)
    obs = env_reset!(env; density_sched = (_ -> DENSITY), year = year)
    done = false; acts = Vector{Vector{Float64}}()
    while !done
        a = vec(Float64.(actor_sample(actor, Float32.(reshape(obs, :, 1)); deterministic = true)[1]))
        push!(acts, a); obs, _, done, _ = env_step!(env, a)
    end
    return reduce(hcat, acts)
end

function block_init(A, ndays)
    M = Matrix{Float64}(undef, ACT_DIM, NB)
    for b in 1:NB
        cols = [d+1 for d in 0:ndays-1 if day_block(d, ndays, NB) == b]
        M[:, b] = isempty(cols) ? A[:, end] : vec(mean(A[:, cols], dims = 2))
    end
    return vec(M)
end

function cem_year(env, year, z0, ndays)   # returns (best_value, best_schedule)
    dim = length(z0); μ = copy(z0); σ = fill(STD0, dim)
    best_z = copy(z0); best_v = eval_block(env, z0, year, ndays)
    for _ in 1:ITERS
        Z = [clamp.(μ .+ σ .* randn(dim), 0.0, 1.0) for _ in 1:POP]
        V = [eval_block(env, z, year, ndays) for z in Z]
        idx = partialsortperm(V, 1:NEL, rev = true)
        E = reduce(hcat, Z[idx]); μ = vec(mean(E, dims = 2)); σ = vec(std(E, dims = 2)) .+ 1e-3
        if V[idx[1]] > best_v; best_v = V[idx[1]]; best_z = Z[idx[1]]; end
    end
    return best_v, best_z
end

function main()
    actor = deserialize(joinpath(@__DIR__, "sac_actor_best.jls"))
    env   = TwinEnv(params, gp, bank; forecast_skill = 1.0, rng = MersenneTwister(1), crop = :cohort)  # realistic forecast in obs
    ndays = env.ndays
    OBS = Vector{Float32}[]; ACT = Vector{Float32}[]
    REW = Float32[]; OB2 = Vector{Float32}[]; DON = Float32[]
    for (k, y) in enumerate(TRAIN_YEARS)
        z0 = block_init(sac_trace(env, actor, y), ndays)
        v, zb = cem_year(env, y, z0, ndays)
        S = reshape(zb, ACT_DIM, NB)
        obs = env_reset!(env; density_sched = (_ -> DENSITY), year = y); done = false; day = 0
        while !done
            a = clamp.(S[:, day_block(day, ndays, NB)], 0.0, 1.0)
            obs2, r, done, _ = env_step!(env, a)
            push!(OBS, Float32.(obs)); push!(ACT, Float32.(a))
            push!(REW, Float32(r)); push!(OB2, Float32.(obs2)); push!(DON, Float32(done))
            obs = obs2; day += 1
        end
        @printf("  year %d (%d/%d): expert profit %.2f\n", y, k, length(TRAIN_YEARS), v); flush(stdout)
    end
    Om = reduce(hcat, OBS); Am = reduce(hcat, ACT); O2 = reduce(hcat, OB2)
    fn = SMOKE ? "mpc_demos_smoke.jls" : "mpc_demos.jls"
    serialize(joinpath(@__DIR__, fn), (obs = Om, act = Am, rew = REW, obs2 = O2, done = DON))
    @printf("saved %d expert transitions -> %s (obs %dx%d, act %dx%d, +rew/obs2/done)\n",
            size(Om, 2), fn, size(Om, 1), size(Om, 2), size(Am, 1), size(Am, 2))
end

main()
