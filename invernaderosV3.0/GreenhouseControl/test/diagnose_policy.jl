# =============================================================================
# diagnose_policy.jl -- does the trained SAC policy TIME lamps to the weather?
#   Loads the saved actor and logs, per day, the lamp decision vs the solar light.
# =============================================================================
using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Flux, Serialization, Dates, Random, Printf, Statistics

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
const DENSITY = 2.90

actor = deserialize(joinpath(@__DIR__, "sac_actor_best.jls"))   # campaign best-VAL winner
"Deterministic policy: normalized action in [0,1]^10 from the actor mean."
det(obs) = begin
    o  = actor(Float32.(reshape(obs, :, 1)))
    mu = o[1:11, 1]
    Float64.((tanh.(mu) .+ 1) ./ 2)
end
daily_solar(env, d) = mean(env.weather.Iglobal(env.start_unix + d*86400.0 + h*3600.0) for h in 0:23)

env = TwinEnv(params, gp, bank; crop=:cohort, forecast_skill=0.0)

# CEM constant reference (what SAC must "beat" by reacting)
CEM_LAMP_H = 3.19

for y in (1990, 2003, 2014, 2015)   # held-out TEST years
    obs = env_reset!(env; density_sched=(_->DENSITY), year=y)
    day=Int[]; solar=Float64[]; lamph=Float64[]; lstart=Float64[]; Tday=Float64[]; co2=Float64[]; prof=Float64[]
    done=false
    while !done
        d = env.day
        a = det(obs)
        phys = ACT_LO .+ a .* (ACT_HI .- ACT_LO)
        s = daily_solar(env, d)
        obs, r, done, _ = env_step!(env, a)
        push!(day,d); push!(solar,s); push!(lamph,phys[4]); push!(lstart,phys[10])
        push!(Tday,phys[1]); push!(co2,phys[3]); push!(prof,r)
    end
    ρ = std(lamph) > 1e-6 ? cor(lamph, solar) : NaN
    @printf("\n===== year %d =====  season profit %.2f EUR/m2\n", y, sum(prof))
    @printf("mean lamp hours %.2f (CEM constant %.2f) | corr(lamp_hours, solar) = %+.2f\n",
            mean(lamph), CEM_LAMP_H, ρ)
    @printf("lamp hours: bright-third days %.2f  vs  dark-third days %.2f\n",
            mean(lamph[sortperm(solar,rev=true)[1:33]]), mean(lamph[sortperm(solar)[1:33]]))
    @printf("%-4s %8s %7s %7s %6s %7s\n","day","solar","lamp_h","l_start","Tday","CO2")
    for i in 1:10:length(day)
        @printf("%-4d %8.1f %7.2f %7.2f %6.1f %7.0f\n", day[i],solar[i],lamph[i],lstart[i],Tday[i],co2[i])
    end
end
