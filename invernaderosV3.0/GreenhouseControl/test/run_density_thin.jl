# =============================================================================
# run_density_thin.jl -- efficient strategy: dense early, thin down; optimize thin day
# =============================================================================
# Run:  julia --project=. test/run_density_thin.jl

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates, Printf

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params = update_params(load_params(climate_json, crop_json), ["psi2","J_max"], [27800.0, 1.15e-4])
gp = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0])

start_date = DateTime(1998, 7, 11, 0, 0); ndays = 100
weather = load_weather(meteo, start_date, start_date + Day(ndays + 5); sky_from_clouds = true)
sp = daynight_setpoints(Tset_day=22+273.15, Tset_night=19+273.15, CO2_set=1200.0,
                        light_start=4.0, light_end=20.0, VentpBand=4.0, ofset=1.0, ToutMax=12+273.15)

function run(tag, ds)
    res = simulate_twin(params, gp, weather, sp, datetime2unix(start_date), ndays; LAI0=0.5, density_sched=ds)
    e = season_economics(res; cropdays=ndays)
    @printf("%-26s yield %5.1f kg  LAImax %.1f  leaf_end %4.0f g  NET %6.1f EUR/m2\n",
            tag, res.yield_FW[end], maximum(res.LAI), res.W_leaf[end], e.net)
end

println("=== references ===")
run("constant 2.5",           density_schedule([(0,2.5)]))
run("constant 4.0",           density_schedule([(0,4.0)]))
println("=== dense 4.0 -> thin to 2.5 at day: (sweep the thinning day) ===")
for td in (8, 12, 16, 20, 25, 30, 40, 55)
    run("dense4.0->2.5 @day$td", density_schedule([(0,4.0),(td,2.5)]))
end
println("=== deeper thin 4.0 -> 2.0 ===")
for td in (16, 25, 35)
    run("dense4.0->2.0 @day$td", density_schedule([(0,4.0),(td,2.0)]))
end
