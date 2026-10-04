# =============================================================================
# diagnose_realism.jl -- fair matched-season check + floor heating isolation
# =============================================================================
# Run:  julia --project=. test/diagnose_realism.jl
#
# Runs the twin on the SAME 2018 autumn season the AGC Reference grower farmed
# (Aug 14, 116 days), so our resource use and yield are comparable to the known
# reference numbers, and contrasts floor-on vs floor-off to see how much the
# thermal-mass anchor adds to heating.

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

gp = load_growth_params(growth_json)
gp = update_params(gp, ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max"],
                       [0.091, 21.0, 2.5, 16.0, 1.5, 16.0])

# 1998 summer season (the long meteo file has gaps: 2018 Aug-Dec is missing).
# This isolates the FLOOR heating contribution and night-setpoint sensitivity.
start_date = DateTime(1998, 7, 11, 0, 0); ndays = 100
weather    = load_weather(meteo, start_date, start_date + Day(ndays + 5); sky_from_clouds = true)
sp = daynight_setpoints(Tset_day = 22+273.15, Tset_night = 19+273.15, CO2_set = 1200.0,
                        light_start = 4.0, light_end = 20.0, VentpBand = 4.0, ofset = 1.0,
                        ToutMax = 12+273.15)

function report(tag, params)
    res = simulate_twin(params, gp, weather, sp, datetime2unix(start_date), ndays; LAI0 = 0.5)
    e   = season_economics(res; cropdays = ndays)
    println("\n[$tag]  yield=", round(res.yield_FW[end],digits=1), " kg/m2",
            "  heat=", round(sum(res.heat_kWh),digits=0), " kWh/m2",
            "  CO2=", round(sum(res.co2_kg),digits=1), " kg/m2",
            "  elec=", round(sum(res.lamp_kWh_peak)+sum(res.lamp_kWh_off),digits=0), " kWh/m2",
            "  NET=", round(e.net,digits=1), " EUR/m2")
    return res
end

base = load_params(climate_json, crop_json)
base = update_params(base, ["psi2","J_max"], [27800.0, 1.15e-4])

println("=== 1998 summer season (Jul 11, 100 d) -- floor & setpoint heating diagnosis ===")
report("twin, floor ON  (k_ground=4.33)", base)
report("twin, floor OFF (k_ground=0)   ", update_params(base, ["k_ground"], [0.0]))

# night-setpoint sensitivity: lower the night target to 17C, floor ON
sp = daynight_setpoints(Tset_day = 22+273.15, Tset_night = 17+273.15, CO2_set = 1200.0,
                        light_start = 4.0, light_end = 20.0, VentpBand = 4.0, ofset = 1.0,
                        ToutMax = 12+273.15)
report("twin, night 17C, floor ON      ", base)
