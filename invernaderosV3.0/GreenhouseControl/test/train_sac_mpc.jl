# =============================================================================
# train_sac_mpc.jl -- MPC->RL distillation campaign. Same pipeline as
# train_sac_long.jl, but the demo buffer + behavior-cloning tether use the
# forecast-driven EXPERT demonstrations (mpc_demos.jls) with STATE-CONDITIONED
# imitation, so the agent learns anticipatory (forecast-using) control instead of
# being tethered to a forecast-blind constant program.
#   smoke:  julia -t auto --project=. test/train_sac_mpc.jl smoke   (needs mpc_demos_smoke.jls)
#   full :  julia -t auto --project=. test/train_sac_mpc.jl         (needs mpc_demos.jls)
# =============================================================================
include(joinpath(@__DIR__, "train_sac.jl"))
using Printf, Random, Statistics, Serialization, Dates

const SMOKE = "smoke" in ARGS

# ---- expert demonstrations (obs, act, rew, obs2, done) ----------------------
DEMO_FILE = joinpath(@__DIR__, SMOKE ? "mpc_demos_smoke.jls" : "mpc_demos.jls")
D = deserialize(DEMO_FILE)
const EXP_OBS  = D.obs; const EXP_ACT = D.act; const EXP_REW = D.rew
const EXP_OBS2 = D.obs2; const EXP_DON = D.done
const EXP_U = Float32.(atanh.(clamp.(2 .* EXP_ACT .- 1f0, -0.999f0, 0.999f0)))  # pre-tanh targets
const NEXP = size(EXP_OBS, 2)
@printf("loaded %d expert transitions from %s\n", NEXP, basename(DEMO_FILE))

# ---- data splits ------------------------------------------------------------
ALL   = sort(bank.years)
TEST  = intersect([1990,1997,2003,2009,2014,2015], ALL)
VAL   = intersect([1992,1999,2005,2011,2016,2017], ALL)
TRAIN = setdiff(ALL, vcat(TEST, VAL))
train_bank = WeatherBank(bank.df, bank.plant_month, bank.plant_day, bank.ndays, TRAIN, bank.sky_from_clouds)
@printf("splits: TRAIN %d, VAL %d, TEST %d years\n", length(TRAIN), length(VAL), length(TEST))

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
function profits(actfn, evenv, years)
    map(years) do y
        obs = env_reset!(evenv; density_sched=(_->DENSITY), year=y); tot=0.0; done=false
        while !done; obs, r, done, _ = env_step!(evenv, actfn(obs)); tot += r; end
        tot
    end
end
det_actfn(actor) = obs -> Float64.((tanh.(actor(Float32.(reshape(obs,:,1)))[1:ACT_DIM,1]) .+ 1) ./ 2)
eval_const(A, evenv, years) = mean(map(years) do y
    obs = env_reset!(evenv; density_sched=(_->DENSITY), year=y); tot=0.0; done=false
    while !done; obs, r, done, _ = env_step!(evenv, A); tot += r; end
    tot
end)

const DEMO_FRAC = 0.25
const CRITIC_PRETRAIN = 4000
const LAMBDA_BC       = 3.0f0
const BC_DECAY        = 40000f0
function mix_minibatch(buf, demo, n::Int, frac::Float64)
    nd = round(Int, frac * n)
    (nd == 0 || nsamp(demo) == 0) && return minibatch(buf, n)
    s1,a1,r1,t1,d1 = minibatch(buf,  n - nd)
    s2,a2,r2,t2,d2 = minibatch(demo, nd)
    return (hcat(s1,s2), hcat(a1,a2), vcat(r1,r2), hcat(t1,t2), vcat(d1,d2))
end

function train_one(seed::Int, cfg::SAC)
    Random.seed!(seed)
    env  = TwinEnv(params, gp, train_bank; forecast_skill=1.0, rng=MersenneTwister(seed),     crop=:cohort)
    evenv= TwinEnv(params, gp, bank;       forecast_skill=0.0, rng=MersenneTwister(seed+100), crop=:cohort)
    ag = build(cfg); buf = Buffer(cfg.buffer_cap); demo = Buffer(30000)

    # demo buffer = protected EXPERT transitions; also seed the online buffer with them
    for j in 1:NEXP
        r = Float32(EXP_REW[j]*cfg.reward_scale)
        store!(demo, EXP_OBS[:,j], EXP_ACT[:,j], r, EXP_OBS2[:,j], EXP_DON[j])
        store!(buf,  EXP_OBS[:,j], EXP_ACT[:,j], r, EXP_OBS2[:,j], EXP_DON[j])
    end
    # BC warm-start: fit actor to EXPERT actions at expert states (state-conditioned)
    for _ in 1:500
        idx = rand(1:NEXP, 128); S = EXP_OBS[:, idx]; U = EXP_U[:, idx]
        g = Flux.gradient(m -> mean(sum((m(S)[1:ACT_DIM,:] .- U).^2, dims=1)), ag.actor)[1]
        Flux.update!(ag.oa, ag.actor, g)
    end
    best_val = eval_actor(ag.actor, evenv, VAL); best_actor = deepcopy(ag.actor)
    @printf("   seed %d warm-start(expert) VAL %.2f\n", seed, best_val); flush(stdout)

    for _ in 1:CRITIC_PRETRAIN
        sac_update!(ag, cfg, mix_minibatch(buf, demo, cfg.batch, DEMO_FRAC); update_actor = false)
    end
    @printf("   seed %d critic-warmup done, VAL %.2f\n", seed, eval_actor(ag.actor, evenv, VAL)); flush(stdout)

    obs = env_reset!(env; density_sched=(_->DENSITY))
    for step in 1:cfg.total_steps
        a = if step <= cfg.warmup
            rand(Float32, ACT_DIM)
        else
            av, _ = actor_sample(ag.actor, Float32.(reshape(obs,:,1))); vec(av)
        end
        obs2, r, done, _ = env_step!(env, Float64.(a))
        store!(buf, Float32.(obs), a, Float32(r*cfg.reward_scale), Float32.(obs2), Float32(done))
        obs = obs2
        done && (obs = env_reset!(env; density_sched=(_->DENSITY)))
        if step > cfg.update_after
            bc_w = LAMBDA_BC * max(0f0, 1f0 - Float32(step)/BC_DECAY)
            idx = rand(1:NEXP, cfg.batch); bS = EXP_OBS[:, idx]; bU = EXP_U[:, idx]
            for _ in 1:cfg.updates_per_step
                sac_update!(ag, cfg, mix_minibatch(buf, demo, cfg.batch, DEMO_FRAC);
                            bc_states=bS, bc_actions_u=bU, bc_weight=bc_w)
            end
        end
        if step % cfg.eval_every == 0
            v = eval_actor(ag.actor, evenv, VAL)
            v > best_val && (best_val = v; best_actor = deepcopy(ag.actor))
            @printf("   seed %d step %6d | VAL %.2f (best %.2f) | alpha %.3f\n",
                    seed, step, v, best_val, exp(ag.logα[1])); flush(stdout)
        end
    end
    return (; seed, val=best_val, test=eval_actor(best_actor, evenv, TEST), best_actor)
end

BUDGET   = SMOKE ? 600.0 : 8*3600.0
MAXSEEDS = SMOKE ? 1 : 15
cfg = SMOKE ?
      SAC(lr=1f-4, total_steps=1500,  eval_every=500,  warmup=200, update_after=200, updates_per_step=2, target_entropy=-3f0) :
      SAC(lr=1f-4, total_steps=100000, eval_every=4000, warmup=0,   update_after=100, updates_per_step=2, target_entropy=-3f0)

evenv0 = TwinEnv(params, gp, bank; forecast_skill=0.0, crop=:cohort)
CEM_A  = Float64.(WARM_A)
cem_test = eval_const(CEM_A, evenv0, TEST); cem_val = eval_const(CEM_A, evenv0, VAL)
@printf("CEM constant: VAL %.2f  TEST %.2f  EUR/m2\n", cem_val, cem_test)
@printf("budget %.1f h, up to %d seeds, %d steps/seed, demo_frac %.2f (EXPERT demos)\n\n",
        BUDGET/3600, MAXSEEDS, cfg.total_steps, DEMO_FRAC); flush(stdout)

t0 = time(); res = NamedTuple[]; seed = 200
while true
    global seed += 1
    r = train_one(seed, cfg); push!(res, r)
    @printf("=== seed %d done: VAL %.2f  TEST %.2f  (CEM TEST %.2f)  [%.0f min]\n\n",
            r.seed, r.val, r.test, cem_test, (time()-t0)/60); flush(stdout)
    mkpath(joinpath(@__DIR__, "campaign_actors_mpc"))
    serialize(joinpath(@__DIR__, "campaign_actors_mpc", "seed_$(r.seed).jls"), r.best_actor)
    r.val >= maximum(x.val for x in res) && serialize(joinpath(@__DIR__, "sac_actor_mpc_best.jls"), r.best_actor)
    el = time()-t0; avg = el/length(res)
    (el + avg > BUDGET || length(res) >= MAXSEEDS) && break
end

tests = Float64[r.test for r in res]
println("\n=========================== MPC-DISTILL CAMPAIGN SUMMARY ===========================")
@printf("%d seeds in %.1f h | steps/seed %d\n", length(res), (time()-t0)/3600, cfg.total_steps)
@printf("CEM constant   : VAL %.2f  TEST %.2f\n", cem_val, cem_test)
@printf("SAC-MPC        : TEST %.2f +/- %.2f  (best seed TEST %.2f)\n",
        mean(tests), length(tests)>1 ? std(tests) : 0.0, maximum(tests))
@printf("beats CEM on TEST in %d / %d seeds ; median %.2f worst %.2f\n",
        count(>(cem_test), tests), length(tests), median(tests), minimum(tests))
@printf("per-seed TEST: %s\n", join((@sprintf("%.2f", t) for t in tests), ", "))
best_seed = res[argmax([r.val for r in res])]
cem_yr = profits(obs->CEM_A, evenv0, TEST); sac_yr = profits(det_actfn(best_seed.best_actor), evenv0, TEST)
println("\nper-TEST-year profit (best-VAL seed $(best_seed.seed)):")
for (i,y) in enumerate(TEST); @printf("%-6d CEM %6.2f   SAC-MPC %6.2f\n", y, cem_yr[i], sac_yr[i]); end
@printf("mean   CEM %6.2f   SAC-MPC %6.2f\n", mean(cem_yr), mean(sac_yr))
open(joinpath(@__DIR__, "campaign_manifest_mpc.csv"), "w") do io
    println(io, "seed,val,test"); for r in res; println(io, "$(r.seed),$(r.val),$(r.test)"); end
end
println("\nbest -> test/sac_actor_mpc_best.jls ; pool -> test/campaign_actors_mpc/ ; manifest -> campaign_manifest_mpc.csv")
