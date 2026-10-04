# =============================================================================
# compare_vs_agc.jl -- trained agent vs the 6 real AGC-2018 cucumber growers,
# with an ELECTRICITY-PRICE SWEEP (the agent's trajectory is fixed; only the cost
# accounting changes, so we score everyone at each price). Shows whether the
# agent's profit edge is real control or just a harsh-price artifact.
#
# Real per-grower data: invernaderos2.0/data/<team>/{Production,ResourceCalculations}.csv
# BOTH growers and agent scored on OUR economics (heat 0.09, CO2 0.30, fixed 15/yr).
# Agent runs on a MATCHED crop (Reference calibrated params, density 2.5) and on its
# NATIVE config. Caveat: agent was TRAINED at high elec price, so at cheap power it
# will NOT light more (fixed policy) -- a truly fair cheap-power test needs retraining.
#
#   julia -t auto --project=. test/compare_vs_agc.jl [actor_file]
# =============================================================================
include(joinpath(@__DIR__, "train_sac.jl"))
using DelimitedFiles, Statistics, Printf, Serialization, Random

const TEST_YEARS = [1990, 1997, 2003, 2009, 2014, 2015]
const PRICES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]   # EUR/kWh electricity sweep
const DATA  = normpath(joinpath(@__DIR__, "..", "..", "..", "invernaderos2.0", "data"))
const TPAR  = normpath(joinpath(@__DIR__, "..", "..", "..", "invernaderosV2.2", "GreenhouseSim", "test", "figures", "teams_params.csv"))
const TEAMS = ["AiCU", "Croperators", "DeepGreens", "Reference(Growers)", "Sonoma", "iGrow"]

const PRICE_CUKE      = GreenhouseControl.PRICE_CUKE
const COST_HEAT       = GreenhouseControl.COST_HEAT
const COST_CO2        = GreenhouseControl.COST_CO2
const FIXED_GH_PER_YR = GreenhouseControl.FIXED_GH_PER_YR
const PLANT_COST_STEM = GreenhouseControl.PLANT_COST_STEM

# net profit at electricity price ep [EUR/kWh]
netprofit(t, dens, ep) = t.yield_kg*PRICE_CUKE - t.heat*COST_HEAT - t.elec*ep -
                         t.co2*COST_CO2 - (t.days/365)*FIXED_GH_PER_YR - PLANT_COST_STEM*dens

function team_totals(team)
    P = readdlm(joinpath(DATA, team, "Production.csv"), ',', skipstart=1)
    R = readdlm(joinpath(DATA, team, "ResourceCalculations.csv"), ',', skipstart=1)
    (; yield_kg = P[end, 6], co2 = sum(R[:,1]), elec = sum(R[:,2]), heat = sum(R[:,3]), days = size(R,1))
end

function load_team_params()
    M = readdlm(TPAR, ',', skipstart=1); d = Dict{String,NamedTuple}()
    for i in 1:size(M,1)
        d[String(M[i,1])] = (; J_max=M[i,2], node_rate=M[i,4], veg_sink_max=M[i,6],
                               LAI_max=M[i,8], set_rate=M[i,10], Wf_max=M[i,12], set_start_day=M[i,14])
    end
    return d
end

function agent_totals(actor, p, g, dens, years)
    env = TwinEnv(p, g, bank; forecast_skill=0.0, rng=MersenneTwister(1), crop=:cohort)
    Y=Float64[]; H=Float64[]; E=Float64[]; C=Float64[]
    for y in years
        obs = env_reset!(env; density_sched=(_->dens), year=y); done=false
        h=0.0; e=0.0; c=0.0
        while !done
            a = vec(Float64.(actor_sample(actor, Float32.(reshape(obs,:,1)); deterministic=true)[1]))
            obs, _, done, info = env_step!(env, a); h += info.hkWh; e += info.elec; c += info.ckg
        end
        push!(Y, env.gs.yield_FW); push!(H,h); push!(E,e); push!(C,c)
    end
    (; yield_kg=mean(Y), heat=mean(H), elec=mean(E), co2=mean(C), days=env.ndays)
end

function main()
    afile = isempty(ARGS) ? "sac_actor_mpc_best.jls" : ARGS[1]
    isfile(joinpath(@__DIR__, afile)) || (afile = "sac_actor_best.jls")
    actor = deserialize(joinpath(@__DIR__, afile)); tp = load_team_params()
    @printf("actor: %s\n\n", afile)

    reals = [(t, team_totals(t), 2.5) for t in TEAMS]
    ref = tp["Reference(Growers)"]
    p_ref = update_params(params, ["J_max"], [ref.J_max])
    g_ref = update_params(gp, ["node_rate","veg_sink_max","LAI_max","set_rate","Wf_max","set_start_day"],
                          [ref.node_rate, ref.veg_sink_max, ref.LAI_max, ref.set_rate, ref.Wf_max, ref.set_start_day])
    rows = vcat(reals,
                [("AGENT ref-crop", agent_totals(actor, p_ref, g_ref, 2.5, TEST_YEARS), 2.5),
                 ("AGENT native",   agent_totals(actor, params, gp, DENSITY, TEST_YEARS), DENSITY)])

    println("resources (season totals; agent = mean over ", length(TEST_YEARS), " autumn years):")
    @printf("%-16s %7s %8s %8s %7s\n", "who", "kg/m2", "heat", "elec", "CO2kg")
    for (n,t,_) in rows; @printf("%-16s %7.1f %8.1f %8.1f %7.1f\n", n, t.yield_kg, t.heat, t.elec, t.co2); end

    println("\nnet profit [EUR/m2] vs electricity price [EUR/kWh]  (heat 0.09, CO2 0.30 fixed):")
    @printf("%-16s", "who"); for p in PRICES; @printf("%8.2f", p); end; println()
    for (n,t,dens) in rows
        @printf("%-16s", n); for p in PRICES; @printf("%8.2f", netprofit(t, dens, p)); end; println()
    end
    println("\nagent policy is FIXED (trained at ~0.25); at cheap power it won't light more.")
    println("A fair cheap-power test would RETRAIN the agent at that price (it would then light + grow more).")
end

main()
