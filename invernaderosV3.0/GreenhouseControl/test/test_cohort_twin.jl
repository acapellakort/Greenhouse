# =============================================================================
# test_cohort_twin.jl -- cohort vs big-leaf crop in the RL env: yield & profit
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
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max","leaf_per_node","leaf_lifespan_dd"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0, 0.040, 600.0])
bank = WeatherBank(meteo; plant_month=8, plant_day=14, ndays=100, sky_from_clouds=true)

norm(phys) = (phys .- ACT_LO) ./ (ACT_HI .- ACT_LO)
COLD = norm([16.28, 16.29, 789.83, 0.0, 6.21, 1.66, 11.53, 0.20, 49.51, 4.0, 1.0])   # CEM optimum (big-leaf)
WARM = norm([23.0,  20.0,  1000.0, 0.0, 4.0,  1.0,  20.0,  0.60, 50.0, 4.0, 1.0])     # warmer, more CO2, no lamps
LIT  = norm([22.0,  20.0,  1000.0, 18.0, 4.0, 1.0,  20.0,  0.60, 50.0, 4.0, 1.0])     # warm + 18h supplemental lamps

function run_season(env, a, year; dens=2.5)
    env_reset!(env; density_sched=(_->dens), year=year)
    tot=0.0; done=false
    while !done; _, r, done, _ = env_step!(env, a); tot += r; end
    return tot, env.gs.yield_FW
end

for crop in (:bigleaf, :cohort)
    env = TwinEnv(params, gp, bank; crop=crop, forecast_skill=0.0)
    println("\n=== crop = $crop ===")
    @printf("%-6s %-6s  %8s  %8s\n", "prog", "year", "yield_kg", "profit")
    for (tag, a) in (("COLD",COLD), ("WARM",WARM), ("LIT",LIT)), y in (1995, 2007)
        p, yld = run_season(env, a, y)
        @printf("%-6s %-6d  %8.1f  %8.2f\n", tag, y, yld, p)
    end
end
