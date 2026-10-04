# =============================================================================
# train_sac_long.jl -- ~6h SAC campaign with TRAIN/VAL/TEST splits, validation
#   best-checkpointing, stability-tuned hyperparameters, and a wall-clock loop.
#   Run:  julia -t auto --project=. test/train_sac_long.jl
# =============================================================================
# Reuses the SAC machinery (Buffer/build/actor_sample/sac_update!/WARM_A/params/
# gp/bank/DENSITY) from train_sac.jl (guarded: including it does NOT auto-train).
include(joinpath(@__DIR__, "train_sac.jl"))
using Printf, Random, Statistics, Serialization, Dates

# ---- data splits (disjoint; the 4 CEM-tuning years stay in TRAIN) -----------
ALL   = sort(bank.years)
TEST  = intersect([1990,1997,2003,2009,2014,2015], ALL)   # never used for training/selection
VAL   = intersect([1992,1999,2005,2011,2016,2017], ALL)   # best-checkpoint selection
TRAIN = setdiff(ALL, vcat(TEST, VAL))
@assert isempty(intersect(TEST, VAL)) && all(y-> y in TRAIN, [1995,2001,2007,2013])
train_bank = WeatherBank(bank.df, bank.plant_month, bank.plant_day, bank.ndays, TRAIN, bank.sky_from_clouds)
@printf("splits: TRAIN %d, VAL %d, TEST %d years\n", length(TRAIN), length(VAL), length(TEST))

# ---- deterministic eval of an actor (Chain) on a year list ------------------
function eval_actor(actor, evenv, years)
    mean(map(years) do y
        obs = env_reset!(evenv; density_sched=(_->DENSITY), year=y); tot=0.0; done=false
        while !done
            o = actor(Float32.(reshape(obs,:,1))); a = Float64.((tanh.(o[1:ACT_DIM,1]) .+ 1) ./ 2)
            obs, r, done, _ = env_step!(evenv, a); tot += r
        end
        tot
    end)
end
# per-year profit vector for a policy function obs->action (for TEST breakdowns)
function profits(actfn, evenv, years)
    map(years) do y
        obs = env_reset!(evenv; density_sched=(_->DENSITY), year=y); tot=0.0; done=false
        while !done; obs, r, done, _ = env_step!(evenv, actfn(obs)); tot += r; end
        tot
    end
end
det_actfn(actor) = obs -> Float64.((tanh.(actor(Float32.(reshape(obs,:,1)))[1:ACT_DIM,1]) .+ 1) ./ 2)

function eval_const(A, evenv, years)
    mean(map(years) do y
        obs = env_reset!(evenv; density_sched=(_->DENSITY), year=y); tot=0.0; done=false
        while !done; obs, r, done, _ = env_step!(evenv, A); tot += r; end
        tot
    end)
end

# ---- demonstration-anchored replay: mix a fraction of each minibatch from a
#      protected CEM/expert buffer (0 = plain replay; ~0.25 = demo-anchored) --------
const DEMO_FRAC = 0.25
const CRITIC_PRETRAIN = 4000     # critic-only updates on the demo buffer before the actor moves
const LAMBDA_BC       = 3.0f0    # initial behavior-cloning tether weight toward the CEM program
const BC_DECAY        = 40000f0  # steps over which the BC tether decays linearly to 0
function mix_minibatch(buf, demo, n::Int, frac::Float64)
    nd = round(Int, frac * n)
    (nd == 0 || nsamp(demo) == 0) && return minibatch(buf, n)
    s1,a1,r1,t1,d1 = minibatch(buf,  n - nd)
    s2,a2,r2,t2,d2 = minibatch(demo, nd)
    return (hcat(s1,s2), hcat(a1,a2), vcat(r1,r2), hcat(t1,t2), vcat(d1,d2))
end

# ---- one seed: warm-start -> train -> keep best-VAL actor -> TEST it ---------
function train_one(seed::Int, cfg::SAC)
    Random.seed!(seed)
    env  = TwinEnv(params, gp, train_bank; forecast_skill=1.0, rng=MersenneTwister(seed),     crop=:cohort)
    evenv= TwinEnv(params, gp, bank;       forecast_skill=0.0, rng=MersenneTwister(seed+100), crop=:cohort)
    ag = build(cfg); buf = Buffer(cfg.buffer_cap)
    demo = Buffer(30000)                                   # protected CEM/expert transitions
    for _ in 1:10
        obs = env_reset!(env; density_sched=(_->DENSITY)); done=false
        while !done
            a = clamp.(WARM_A .+ 0.03f0 .* randn(Float32, ACT_DIM), 0f0, 1f0)
            obs2, r, done, _ = env_step!(env, Float64.(a))
            store!(demo, Float32.(obs), a, Float32(r*cfg.reward_scale), Float32.(obs2), Float32(done))
            obs = obs2
        end
    end

    # warm start: seed buffer with CEM rollouts + behavior-clone actor to CEM program
    bc = Vector{Float32}[]
    for _ in 1:6
        obs = env_reset!(env; density_sched=(_->DENSITY)); done=false
        while !done
            a = clamp.(WARM_A .+ 0.03f0 .* randn(Float32, ACT_DIM), 0f0, 1f0)
            obs2, r, done, _ = env_step!(env, Float64.(a))
            store!(buf, Float32.(obs), a, Float32(r*cfg.reward_scale), Float32.(obs2), Float32(done))
            push!(bc, Float32.(obs)); obs = obs2
        end
    end
    u_t = Float32.(atanh.(clamp.(2 .* WARM_A .- 1f0, -0.999f0, 0.999f0)))
    for _ in 1:500
        S = reduce(hcat, bc[rand(1:length(bc),128)])
        g = Flux.gradient(m -> mean(sum((m(S)[1:ACT_DIM,:] .- u_t).^2, dims=1)), ag.actor)[1]
        Flux.update!(ag.oa, ag.actor, g)
    end

    best_val = eval_actor(ag.actor, evenv, VAL); best_actor = deepcopy(ag.actor)
    @printf("   seed %d warm-start VAL %.2f\n", seed, best_val); flush(stdout)

    # critic-only warmup: fit Q around the behavior-cloned CEM policy before the actor moves,
    # so the actor improves FROM CEM instead of bolting against an untrained critic
    for _ in 1:CRITIC_PRETRAIN
        sac_update!(ag, cfg, mix_minibatch(buf, demo, cfg.batch, DEMO_FRAC); update_actor = false)
    end
    @printf("   seed %d critic-warmup done (%d steps), VAL %.2f\n", seed, CRITIC_PRETRAIN,
            eval_actor(ag.actor, evenv, VAL)); flush(stdout)

    obs = env_reset!(env; density_sched=(_->DENSITY))
    for step in 1:cfg.total_steps
        a = if step <= cfg.warmup
            rand(Float32, ACT_DIM)
        else
            av, _ = actor_sample(ag.actor, Float32.(reshape(obs,:,1)))
            vec(av)
        end
        obs2, r, done, _ = env_step!(env, Float64.(a))
        store!(buf, Float32.(obs), a, Float32(r*cfg.reward_scale), Float32.(obs2), Float32(done))
        obs = obs2
        done && (obs = env_reset!(env; density_sched=(_->DENSITY)))
        if step > cfg.update_after
            bc_w = LAMBDA_BC * max(0f0, 1f0 - Float32(step) / BC_DECAY)   # tether decays to 0 by BC_DECAY
            for _ in 1:cfg.updates_per_step
                sac_update!(ag, cfg, mix_minibatch(buf, demo, cfg.batch, DEMO_FRAC); bc_target = u_t, bc_weight = bc_w)
            end
        end
        if step % cfg.eval_every == 0
            v = eval_actor(ag.actor, evenv, VAL)
            if v > best_val; best_val = v; best_actor = deepcopy(ag.actor); end
            @printf("   seed %d step %6d | VAL %.2f (best %.2f) | alpha %.3f\n",
                    seed, step, v, best_val, exp(ag.logα[1])); flush(stdout)
        end
    end
    test = eval_actor(best_actor, evenv, TEST)
    return (; seed, val=best_val, test, best_actor)
end

# ---- wall-clock campaign ----------------------------------------------------
# smoke test:  julia -t auto --project=. test/train_sac_long.jl smoke   (~2 min, 1 seed)
const SMOKE    = "smoke" in ARGS
BUDGET   = SMOKE ? 600.0 : 8*3600.0
MAXSEEDS = SMOKE ? 1 : 15
cfg = SMOKE ?
      SAC(lr=1f-4, total_steps=1500,  eval_every=500,  warmup=200, update_after=200, updates_per_step=2, target_entropy=-3f0) :
      SAC(lr=1f-4, total_steps=100000, eval_every=4000, warmup=0,   update_after=100, updates_per_step=2, target_entropy=-3f0)

evenv0 = TwinEnv(params, gp, bank; forecast_skill=0.0, crop=:cohort)
CEM_A  = Float64.(WARM_A)
cem_test = eval_const(CEM_A, evenv0, TEST); cem_val = eval_const(CEM_A, evenv0, VAL)
@printf("CEM constant: VAL %.2f  TEST %.2f  EUR/m2\n", cem_val, cem_test)
@printf("budget %.1f h, up to %d seeds, %d steps/seed, lr %.0e, target-entropy %.1f, demo_frac %.2f\n\n",
        BUDGET/3600, MAXSEEDS, cfg.total_steps, cfg.lr, cfg.target_entropy, DEMO_FRAC); flush(stdout)

t0 = time(); res = NamedTuple[]; seed = 100
while true
    global seed += 1
    r = train_one(seed, cfg); push!(res, r)
    @printf("=== seed %d done: VAL %.2f  TEST %.2f  (CEM TEST %.2f)  [%.0f min elapsed]\n\n",
            r.seed, r.val, r.test, cem_test, (time()-t0)/60); flush(stdout)
    mkpath(joinpath(@__DIR__, "campaign_actors"))
    serialize(joinpath(@__DIR__, "campaign_actors", "seed_$(r.seed).jls"), r.best_actor)  # pool for ensembling
    if r.val >= maximum(x.val for x in res)
        serialize(joinpath(@__DIR__, "sac_actor_best.jls"), r.best_actor)
    end
    el = time()-t0; avg = el/length(res)
    (el + avg > BUDGET || length(res) >= MAXSEEDS) && break
end

# ---- summary ----------------------------------------------------------------
tests = Float64[r.test for r in res]; vals = Float64[r.val for r in res]
println("\n=========================== CAMPAIGN SUMMARY ===========================")
@printf("%d seeds in %.1f h | steps/seed %d\n", length(res), (time()-t0)/3600, cfg.total_steps)
@printf("CEM constant   : VAL %.2f  TEST %.2f\n", cem_val, cem_test)
@printf("SAC best-ckpt  : TEST %.2f +/- %.2f  (best seed TEST %.2f)\n",
        mean(tests), length(tests)>1 ? std(tests) : 0.0, maximum(tests))
@printf("SAC beats CEM on TEST in %d / %d seeds\n", count(>(cem_test), tests), length(tests))
@printf("SAC TEST  median %.2f  worst-seed %.2f  best-seed %.2f\n", median(tests), minimum(tests), maximum(tests))
@printf("per-seed TEST: %s\n", join((@sprintf("%.2f", t) for t in tests), ", "))

# per-TEST-year breakdown: does the best policy rescue the catastrophic year?
best_seed = res[argmax([r.val for r in res])]
cem_yr = profits(obs->CEM_A, evenv0, TEST)
sac_yr = profits(det_actfn(best_seed.best_actor), evenv0, TEST)
println("\nper-TEST-year profit (best-VAL seed $(best_seed.seed)):")
@printf("%-6s %8s %8s\n", "year", "CEM", "SAC")
for (i,y) in enumerate(TEST); @printf("%-6d %8.2f %8.2f\n", y, cem_yr[i], sac_yr[i]); end
@printf("%-6s %8.2f %8.2f   (mean)\n", "", mean(cem_yr), mean(sac_yr))
open(joinpath(@__DIR__, "campaign_manifest.csv"), "w") do io
    println(io, "seed,val,test")
    for r in res; println(io, "$(r.seed),$(r.val),$(r.test)"); end
end
println("\nbest actor -> test/sac_actor_best.jls ; per-seed actors -> test/campaign_actors/ ; manifest -> test/campaign_manifest.csv")
