# =============================================================================
# train_sac.jl -- Soft Actor-Critic (hand-rolled, Flux) on the greenhouse twin
# =============================================================================
# Requires Flux in the project:  julia --project=. -e 'using Pkg; Pkg.add("Flux")'
# Then:                          julia -t auto --project=. test/train_sac.jl
#
# SAC = twin Q-critics + squashed-Gaussian actor + auto entropy temperature.
# Action = 9 normalized setpoints (density fixed here; can be added later).
# Reward = daily net profit EUR/m2. Target to beat: CEM baseline -0.29 EUR/m2.

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Flux, Dates, Random, Printf, Statistics, Serialization

# ---------------------------------------------------------------- environment
INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params = update_params(load_params(climate_json, crop_json), ["psi2","J_max"], [27800.0, 1.15e-4])
gp = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max","leaf_per_node","leaf_lifespan_dd"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0, 0.040, 600.0])
bank = WeatherBank(meteo; plant_month=8, plant_day=14, ndays=100, sky_from_clouds=true)

const OBS_DIM = 18
const ACT_DIM = 11
const DENSITY = 2.40                   # cohort CEM inc1 optimum (lamp+flue-gas)
const WARM_PHYS = [16.01, 15.81, 782.73, 4.00, 4.26, 2.02, 11.75, 0.21, 52.11, 3.98, 0.81]  # CEM inc1 (lamp dimming + flue-gas CO2): mean profit -0.50 EUR/m2
const WARM_A    = Float32.((WARM_PHYS .- ACT_LO) ./ (ACT_HI .- ACT_LO))
const EVAL_YEARS = [1995, 2001, 2007, 2013]

# ---------------------------------------------------------------- SAC config
Base.@kwdef struct SAC
    γ::Float32 = 0.99f0
    τ::Float32 = 0.005f0
    lr::Float32 = 3f-4
    batch::Int = 256
    hidden::Int = 256
    buffer_cap::Int = 100_000
    warmup::Int = 300           # random-action steps (buffer pre-seeded from CEM)
    update_after::Int = 300
    updates_per_step::Int = 2
    total_steps::Int = 24000
    eval_every::Int = 4000
    reward_scale::Float32 = 10f0
    target_entropy::Float32 = -Float32(ACT_DIM)
end

mlp(inp, out, h) = Chain(Dense(inp, h, relu), Dense(h, h, relu), Dense(h, out))

# ---------------------------------------------------------------- replay buffer
mutable struct Buffer
    s::Matrix{Float32}; a::Matrix{Float32}; r::Vector{Float32}
    s2::Matrix{Float32}; d::Vector{Float32}
    cap::Int; idx::Int; full::Bool
end
Buffer(cap) = Buffer(zeros(Float32,OBS_DIM,cap), zeros(Float32,ACT_DIM,cap), zeros(Float32,cap),
                     zeros(Float32,OBS_DIM,cap), zeros(Float32,cap), cap, 0, false)
function store!(b::Buffer, s,a,r,s2,d)
    i = b.idx + 1
    b.s[:,i].=s; b.a[:,i].=a; b.r[i]=r; b.s2[:,i].=s2; b.d[i]=d
    b.idx = i % b.cap; b.full = b.full || i==b.cap
end
nsamp(b::Buffer) = b.full ? b.cap : b.idx
function minibatch(b::Buffer, n)
    N = nsamp(b); ii = rand(1:N, n)
    return b.s[:,ii], b.a[:,ii], b.r[ii], b.s2[:,ii], b.d[ii]
end

# ---------------------------------------------------------------- squashed gaussian
const LOG2 = log(2f0)
"Sample action in [0,1]^9 from the actor + its log-prob (columns = batch)."
function actor_sample(actor, s; deterministic=false)
    o = actor(s)
    μ = o[1:ACT_DIM, :]; logσ = clamp.(o[ACT_DIM+1:end, :], -5f0, 2f0)
    σ = exp.(logσ)
    u = deterministic ? μ : μ .+ σ .* randn(Float32, size(μ))
    t = tanh.(u)                                   # (-1,1)
    a = (t .+ 1f0) .* 0.5f0                          # (0,1)
    # log prob with tanh + affine(/2) correction, summed over action dims
    gl = -0.5f0 .* (((u .- μ) ./ σ).^2) .- logσ .- 0.5f0*log(2f0*π)
    logp = sum(gl .- log.(1f0 .- t.^2 .+ 1f-6) .- LOG2, dims=1)
    return a, logp
end

# ---------------------------------------------------------------- build agent
function build(cfg::SAC)
    actor = mlp(OBS_DIM, 2*ACT_DIM, cfg.hidden)
    q1 = mlp(OBS_DIM+ACT_DIM, 1, cfg.hidden); q2 = mlp(OBS_DIM+ACT_DIM, 1, cfg.hidden)
    q1t = deepcopy(q1); q2t = deepcopy(q2)
    logα = [0f0]
    return (; actor, q1, q2, q1t, q2t, logα,
              oa = Flux.setup(Adam(cfg.lr), actor),
              o1 = Flux.setup(Adam(cfg.lr), q1),
              o2 = Flux.setup(Adam(cfg.lr), q2),
              oα = Flux.setup(Adam(cfg.lr), logα))
end

function soft_update!(targ, src, τ)
    Flux.fmap(targ, src) do pt, ps
        pt isa AbstractArray{<:AbstractFloat} && (pt .= (1-τ).*pt .+ τ.*ps)
        pt
    end
    return nothing
end

function sac_update!(ag, cfg, batch; update_actor::Bool = true, bc_target = nothing, bc_weight::Float32 = 0f0, bc_states = nothing, bc_actions_u = nothing)
    s,a,r,s2,d = batch
    α = exp(ag.logα[1])
    # --- critic target ---
    a2, logp2 = actor_sample(ag.actor, s2)
    q1t = ag.q1t(vcat(s2,a2)); q2t = ag.q2t(vcat(s2,a2))
    minqt = min.(q1t, q2t)
    y = r' .+ cfg.γ .* (1f0 .- d') .* (minqt .- α .* logp2)     # (1,batch)
    # --- critic loss ---
    gq1 = Flux.gradient(m -> mean((m(vcat(s,a)) .- y).^2), ag.q1)[1]
    Flux.update!(ag.o1, ag.q1, gq1)
    gq2 = Flux.gradient(m -> mean((m(vcat(s,a)) .- y).^2), ag.q2)[1]
    Flux.update!(ag.o2, ag.q2, gq2)
    # --- actor loss (+ optional decaying behavior-cloning tether toward CEM) ---
    if update_actor
        ga = Flux.gradient(m -> begin
                ap, logp = actor_sample(m, s)
                q = min.(ag.q1(vcat(s,ap)), ag.q2(vcat(s,ap)))
                base = mean(α .* logp .- q)
                if bc_weight > 0f0
                    if bc_states !== nothing          # state-conditioned imitation of an expert
                        μd = m(bc_states)[1:ACT_DIM, :]
                        base += bc_weight * mean(sum((μd .- bc_actions_u).^2, dims=1))
                    else                               # constant-program tether (legacy)
                        μ = m(s)[1:ACT_DIM, :]
                        base += bc_weight * mean(sum((μ .- bc_target).^2, dims=1))
                    end
                end
                base
            end, ag.actor)[1]
        Flux.update!(ag.oa, ag.actor, ga)
        # --- temperature ---
        _, logp = actor_sample(ag.actor, s)
        gα = Flux.gradient(la -> -mean(la[1] .* (logp .+ cfg.target_entropy)), ag.logα)[1]
        Flux.update!(ag.oα, ag.logα, gα)
    end
    # --- targets ---
    soft_update!(ag.q1t, ag.q1, cfg.τ); soft_update!(ag.q2t, ag.q2, cfg.τ)
    return nothing
end

# ---------------------------------------------------------------- eval
function evaluate(ag, env)
    tot = 0.0
    for y in EVAL_YEARS
        obs = env_reset!(env; density_sched=(_->DENSITY), year=y)
        done=false
        while !done
            a, _ = actor_sample(ag.actor, Float32.(reshape(obs,:,1)); deterministic=true)
            obs, r, done, _ = env_step!(env, vec(Float64.(a)))
            tot += r
        end
    end
    return tot/length(EVAL_YEARS)
end

# ---------------------------------------------------------------- train
function run_seed(seed::Int=2026; verbose::Bool=true, save::Bool=true)
    cfg = SAC()
    Random.seed!(seed)
    env  = TwinEnv(params, gp, bank; forecast_skill=1.0, rng=MersenneTwister(seed), crop=:cohort)
    evenv= TwinEnv(params, gp, bank; forecast_skill=0.0, rng=MersenneTwister(seed+100), crop=:cohort)
    ag   = build(cfg)
    buf  = Buffer(cfg.buffer_cap)

    # ---- WARM START: seed buffer with CEM rollouts + behavior-clone the actor ----
    bc_obs = Vector{Float32}[]
    for _ in 1:6
        obs = env_reset!(env; density_sched=(_->DENSITY)); done=false
        while !done
            a = clamp.(WARM_A .+ 0.03f0 .* randn(Float32, ACT_DIM), 0f0, 1f0)
            obs2, r, done, _ = env_step!(env, Float64.(a))
            store!(buf, Float32.(obs), a, Float32(r*cfg.reward_scale), Float32.(obs2), Float32(done))
            push!(bc_obs, Float32.(obs)); obs = obs2
        end
    end
    u_t = Float32.(atanh.(clamp.(2 .* WARM_A .- 1f0, -0.999f0, 0.999f0)))   # pre-tanh target
    for _ in 1:500
        S = reduce(hcat, bc_obs[rand(1:length(bc_obs), 128)])
        g = Flux.gradient(m -> mean(sum((m(S)[1:ACT_DIM, :] .- u_t).^2, dims=1)), ag.actor)[1]
        Flux.update!(ag.oa, ag.actor, g)
    end
    verbose && @printf("warm-started; pre-RL eval = %.2f EUR/m2 (CEM = -0.77)\n", evaluate(ag, evenv))

    verbose && @printf("SAC[seed %d]: %d steps, warmup %d. cohort CEM baseline -0.77\n", seed, cfg.total_steps, cfg.warmup)
    obs = env_reset!(env; density_sched=(_->DENSITY))
    epret = 0.0; epn = 0
    for step in 1:cfg.total_steps
        a = if step <= cfg.warmup
            rand(Float32, ACT_DIM)
        else
            av,_ = actor_sample(ag.actor, Float32.(reshape(obs,:,1))); vec(av)
        end
        obs2, r, done, _ = env_step!(env, Float64.(a))
        store!(buf, Float32.(obs), a, Float32(r*cfg.reward_scale), Float32.(obs2), Float32(done))
        obs = obs2; epret += r
        if done
            epn += 1
            obs = env_reset!(env; density_sched=(_->DENSITY))
            (verbose && epn % 5 == 0) && @printf("  [ep %3d] train-return %.2f EUR/m2\n", epn, epret)
            epret = 0.0
        end
        if step > cfg.update_after
            for _ in 1:cfg.updates_per_step; sac_update!(ag, cfg, minibatch(buf, cfg.batch)); end
        end
        if step % cfg.eval_every == 0
            ev = evaluate(ag, evenv)
            verbose && @printf("step %5d | eval profit %.2f EUR/m2 | alpha %.3f\n", step, ev, exp(ag.logα[1]))
        end
    end
    ev = evaluate(ag, evenv)
    verbose && @printf("\nFINAL eval profit = %.2f EUR/m2  (cohort CEM inc1 baseline -0.50)\n", ev)

    # ---------------- generalization: CEM fixed program vs SAC, all 31 years -----
    cem_phys = copy(WARM_PHYS)  # CEM -0.50 reference (11-action, inc1)
    CEM_A    = (cem_phys .- ACT_LO) ./ (ACT_HI .- ACT_LO)
    sac_act(obs) = vec(Float64.(actor_sample(ag.actor, Float32.(reshape(obs,:,1)); deterministic=true)[1]))

    function eval_years(actfn, dens, years)
        map(years) do y
            obs = env_reset!(evenv; density_sched=(_->dens), year=y)
            tot=0.0; done=false
            while !done; obs, r, done, _ = env_step!(evenv, actfn(obs)); tot += r; end
            tot
        end
    end

    allyrs = sort(evenv.bank.years)
    tuned  = EVAL_YEARS
    held   = setdiff(allyrs, tuned)
    cem = Dict(zip(allyrs, eval_years(_->CEM_A, DENSITY, allyrs)))
    sac = Dict(zip(allyrs, eval_years(sac_act, DENSITY, allyrs)))

    grp(d, ys) = mean(d[y] for y in ys)
    if verbose
        println("\n=== generalization: mean net profit EUR/m2 (forecast skill=0) ===")
        @printf("%-26s %8s %8s\n", "", "CEM", "SAC")
        @printf("%-26s %8.2f %8.2f\n", "4 CEM-tuned years",  grp(cem,tuned), grp(sac,tuned))
        @printf("%-26s %8.2f %8.2f\n", "27 held-out years",  grp(cem,held),  grp(sac,held))
        @printf("%-26s %8.2f %8.2f\n", "all 31 years",       grp(cem,allyrs),grp(sac,allyrs))
        @printf("%-26s %8.2f %8.2f\n", "worst single year",  minimum(values(cem)), minimum(values(sac)))
    end
    if save
        try
            serialize(joinpath(@__DIR__, "sac_actor.jls"), ag.actor)
            verbose && println("\nsaved trained actor -> test/sac_actor.jls")
        catch e; @warn "actor save failed" e; end
    end
    return (; final_eval = ev,
              sac_tuned=grp(sac,tuned), sac_held=grp(sac,held), sac_all=grp(sac,allyrs), sac_worst=minimum(values(sac)),
              cem_tuned=grp(cem,tuned), cem_held=grp(cem,held), cem_all=grp(cem,allyrs), cem_worst=minimum(values(cem)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_seed()
end
