# =============================================================================
# run_deleaf.jl -- de-leafing: dense stems + shed shaded leaves to a working LAI
# =============================================================================
# Run:  julia --project=. test/run_deleaf.jl

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

function run(tag, dsched, deleaf)
    res = simulate_twin(params, gp, weather, sp, datetime2unix(start_date), ndays;
                        LAI0=0.5, density_sched=dsched, deleaf_sched=deleaf)
    e = season_economics(res; cropdays=ndays)
    @printf("%-34s yield %5.1f kg  LAImax %.1f  leaf_end %4.0f g  NET %6.1f\n",
            tag, res.yield_FW[end], maximum(res.LAI), res.W_leaf[end], e.net)
end

println("=== references (no de-leafing) ===")
run("constant 2.5",              density_schedule([(0,2.5)]), nothing)
run("constant 4.0",              density_schedule([(0,4.0)]), nothing)
println("=== dense stems 4.0 + DE-LEAF to working LAI (keeps stems/fruit) ===")
for D in (12, 20, 30)
    run("dense4.0 + deleaf->2.5 @day$D", density_schedule([(0,4.0)]), density_schedule([(0,99.0),(D,2.5)]))
end
run("dense4.0 + deleaf->3.0 @day20", density_schedule([(0,4.0)]), density_schedule([(0,99.0),(20,3.0)]))
println("=== push it: dense 6.0 stems + deleaf ===")
run("dense6.0 + deleaf->2.5 @day15", density_schedule([(0,6.0)]), density_schedule([(0,99.0),(15,2.5)]))
run("dense6.0 + deleaf->3.0 @day15", density_schedule([(0,6.0)]), density_schedule([(0,99.0),(15,3.0)]))
