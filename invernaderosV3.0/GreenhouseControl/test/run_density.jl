# =============================================================================
# run_density.jl -- variable in-season plant density demo
# =============================================================================
# Run:  julia --project=. test/run_density.jl
# Compares a constant density vs an interplanting schedule (raise density mid-season).

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates

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
                        light_start=4.0, light_end=20.0, VentpBand=4.0, ofset=1.0,
                        ToutMax=12+273.15)

function report(tag, dsched)
    res = simulate_twin(params, gp, weather, sp, datetime2unix(start_date), ndays;
                        LAI0=0.5, density_sched=dsched)
    e = season_economics(res; cropdays=ndays)
    println("[$tag]  density ", round(res.density[1],digits=1), "->", round(res.density[end],digits=1),
            "  LAI ", round(res.LAI[end],digits=2),
            "  yield ", round(res.yield_FW[end],digits=1), " kg/m2",
            "  NET ", round(e.net,digits=1), " EUR/m2")
end

report("constant 2.5      ", density_schedule([(0, 2.5)]))
report("interplant 2.0->3.5", density_schedule([(0, 2.0), (25, 2.75), (55, 3.5)]))
report("dense 3.5          ", density_schedule([(0, 3.5)]))
