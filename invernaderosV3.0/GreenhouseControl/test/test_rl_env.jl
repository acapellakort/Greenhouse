# =============================================================================
# test_rl_env.jl -- smoke test: weather-year sampling + noisy forecast + env
# =============================================================================
using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates, Random, Printf

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params = update_params(load_params(climate_json, crop_json), ["psi2","J_max"], [27800.0, 1.15e-4])
gp = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0])

# ---- weather bank: which years have a fully-covered Aug-14 + 100-day window? --
bank = WeatherBank(meteo; plant_month = 8, plant_day = 14, ndays = 100, sky_from_clouds = true)
@printf("covered start-years (%d): %s\n", length(bank.years), string(bank.years))

env = TwinEnv(params, gp, bank; forecast_skill = 1.0, rng = MersenneTwister(1))

# baseline policy (normalized): Tday22,Tnight19,CO2 1200,lh16,pBand4,ofs1,ToutMax12,VPDmin0.5,Screen50
base_phys = [22.0, 19.0, 1200.0, 16.0, 4.0, 1.0, 12.0, 0.5, 50.0]
a = (base_phys .- ACT_LO) ./ (ACT_HI .- ACT_LO)

obs0 = env_reset!(env; year = bank.years[1])
@printf("\nobs length = %d  (6 crop + 3 climate + %d forecast)\n", length(obs0), 3*3)
println("obs[1..6] crop     = ", round.(obs0[1:6], digits=3))
println("obs[7..9] climate  = ", round.(obs0[7:9], digits=3))
println("obs[10..] forecast = ", round.(obs0[10:end], digits=3))

# ---- run 3 episodes on 3 different sampled years -----------------------------
println("\n=== 3 episodes, fixed baseline policy, domain-randomized weather ===")
@printf("%-6s %8s %8s %8s %8s %8s\n","year","profit","yield","heat","co2","elec")
for ep in 1:3
    env_reset!(env)
    total=0.0; H=0.0; C=0.0; E=0.0; done=false; info=NamedTuple()
    while !done
        _, r, done, info = env_step!(env, a)
        total += r; H += info.hkWh; C += info.ckg; E += info.elec
    end
    @printf("%-6d %8.2f %8.1f %8.1f %8.2f %8.1f\n", info.year, total, env.gs.yield_FW, H, C, E)
end

# ---- forecast skill knob: perfect vs nominal vs useless (same year) ----------
println("\n=== forecast noise sanity (year $(bank.years[1]), lead-3 Tout error) ===")
for sk in (0.0, 1.0, 3.0)
    e2 = TwinEnv(params, gp, bank; forecast_skill = sk, rng = MersenneTwister(7))
    env_reset!(e2; year = bank.years[1])
    fc = forecast_row(e2.forecast, 10)   # day 10 forecast triple-per-lead
    @printf("skill=%.1f  forecast(Tout) leads1-3 = %s\n", sk, string(round.(fc[1:3:end], digits=2)))
end
