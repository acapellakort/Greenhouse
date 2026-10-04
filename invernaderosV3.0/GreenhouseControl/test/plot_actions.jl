# =============================================================================
# plot_actions.jl -- how differently does the agent ACT vs the real Reference
# grower? Rolls out the agent on the MATCHED Reference crop (density 2.5) for one
# autumn year and overlays daily trajectories against the grower's realized data.
# Caveat: agent runs on our historical autumn weather (default 2009), grower on
# 2018 Bleiswijk -- so this shows STRATEGY differences, not a same-weather match.
#
#   julia -t auto --project=. test/plot_actions.jl [actor_file] [year]
#   -> test/figures/agent_vs_grower.png
# =============================================================================
include(joinpath(@__DIR__, "train_sac.jl"))
using DelimitedFiles, Statistics, Printf, Serialization, Random, Plots

const DATA = normpath(joinpath(@__DIR__, "..","..","..","invernaderos2.0","data","Reference(Growers)"))
const TPAR = normpath(joinpath(@__DIR__, "..","..","..","invernaderosV2.2","GreenhouseSim","test","figures","teams_params.csv"))
nm(v) = (w = filter(isfinite, v); isempty(w) ? NaN : mean(w))

function grower_daily()
    G = readdlm(joinpath(DATA,"Greenhouse_climate.csv"), ',', skipstart=1)  # AssimLight1 CO2air3 GHtime5 PipeGrow7 Tair10
    day = floor.(Int, Float64.(G[:,5])) .- 43326
    ds = sort(unique(day)); ds = ds[ds .>= 0]
    Tair = [nm(G[day.==d,10]) for d in ds]
    CO2  = [nm(G[day.==d, 3]) for d in ds]
    lamp = [nm(G[day.==d, 1]) for d in ds]     # AssimLight on/off fraction proxy
    R = readdlm(joinpath(DATA,"ResourceCalculations.csv"), ',', skipstart=1)  # CO2_dosage1 Elec2 Heat3
    P = readdlm(joinpath(DATA,"Production.csv"), ',', skipstart=1)            # Total_Prod_cum6
    (; day=ds, Tair, CO2, lamp, elec=Float64.(R[:,2]), heat=Float64.(R[:,3]), cumyield=Float64.(P[:,6]))
end

function agent_daily(actor, p, g, dens, year)
    env = TwinEnv(p, g, bank; forecast_skill=0.0, rng=MersenneTwister(1), crop=:cohort)
    obs = env_reset!(env; density_sched=(_->dens), year=year); done=false
    T=Float64[]; C=Float64[]; H=Float64[]; E=Float64[]; L=Float64[]; LAI=Float64[]; CY=Float64[]
    Tds=Float64[]; Tns=Float64[]; lint=Float64[]
    span = ACT_HI .- ACT_LO
    while !done
        a = vec(Float64.(actor_sample(actor, Float32.(reshape(obs,:,1)); deterministic=true)[1]))
        phys = ACT_LO .+ a .* span                       # decode to physical setpoints
        obs, _, done, info = env_step!(env, a)
        push!(T, env.last_T); push!(C, env.last_CO2/1.83); push!(H, info.hkWh); push!(E, info.elec)
        push!(LAI, info.LAI); push!(CY, env.gs.yield_FW)
        push!(Tds, phys[1]); push!(Tns, phys[2]); push!(lint, phys[11]); push!(L, info.elec)
    end
    nd = length(T)
    Tout = [mean(env.weather.Tout(env.start_unix + d*86400.0 + h*3600.0) - 273.15 for h in 0:23) for d in 0:nd-1]
    (; day=0:nd-1, Tair=T, Tout, CO2=C, heat=H, elec=E, LAI, cumyield=CY, Tset_day=Tds, Tset_night=Tns, lamp_int=lint)
end

function main()
    afile = isempty(ARGS) ? "sac_actor_mpc_best.jls" : ARGS[1]
    isfile(joinpath(@__DIR__,afile)) || (afile="sac_actor_best.jls")
    year = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2009
    actor = deserialize(joinpath(@__DIR__, afile))
    M = readdlm(TPAR, ',', skipstart=1); ri = findfirst(==("Reference(Growers)"), String.(M[:,1]))
    p_ref = update_params(params, ["J_max"], [M[ri,2]])
    g_ref = update_params(gp, ["node_rate","veg_sink_max","LAI_max","set_rate","Wf_max","set_start_day"],
                          [M[ri,4],M[ri,6],M[ri,8],M[ri,10],M[ri,12],M[ri,14]])
    mode = length(ARGS) >= 3 ? ARGS[3] : "matched"
    A = mode == "native" ? agent_daily(actor, params, gp, DENSITY, year) :
                           agent_daily(actor, p_ref, g_ref, 2.5, year)
    G = grower_daily()
    @printf("crop mode: %s\n", mode)

    kw = (lw=2, legend=:best, grid=true)
    p1 = plot(title="air temperature [C]"; kw...)
        plot!(p1, G.day, G.Tair, label="Reference grower"); plot!(p1, A.day, A.Tair, label="agent")
        plot!(p1, A.day, A.Tset_day, ls=:dash, label="agent Tset_day"); plot!(p1, A.day, A.Tset_night, ls=:dot, label="agent Tset_night")
        plot!(p1, A.day, A.Tout, lw=1, color=:gray, label="outside (agent yr)")
    p2 = plot(title="CO2 [ppm]"; kw...); plot!(p2, G.day, G.CO2, label="grower"); plot!(p2, A.day, A.CO2, label="agent")
    p3 = plot(title="lamp electricity [kWh/m2/day]"; kw...); plot!(p3, G.day, G.elec, label="grower"); plot!(p3, A.day, A.elec, label="agent")
    p4 = plot(title="heating [kWh/m2/day]"; kw...); plot!(p4, G.day, G.heat, label="grower"); plot!(p4, A.day, A.heat, label="agent")
    p5 = plot(title="cumulative yield [kg/m2]"; kw...); plot!(p5, G.day, G.cumyield, label="grower"); plot!(p5, A.day, A.cumyield, label="agent")
    p6 = plot(title="agent lighting program"; kw...); plot!(p6, A.day, A.lamp_int, label="lamp intensity"); plot!(p6, A.day, A.LAI, label="LAI")

    plt = plot(p1,p2,p3,p4,p5,p6, layout=(3,2), size=(1100,1100),
               plot_title="agent ($mode crop, $year weather) vs real Reference grower (2018)")
    mkpath(joinpath(@__DIR__,"figures"))
    out = joinpath(@__DIR__,"figures","agent_vs_grower_$(mode).png"); savefig(plt, out)
    @printf("saved %s\n", out)
    @printf("season totals: agent yield %.1f heat %.0f elec %.0f | grower yield %.1f heat %.0f elec %.0f\n",
            A.cumyield[end], sum(A.heat), sum(A.elec), G.cumyield[end], sum(G.heat), sum(G.elec))
end
main()
