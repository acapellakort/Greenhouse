# =============================================================================
# run_baseline.jl  --  net profit of the current FIXED setpoints (Stage-3 baseline)
# =============================================================================
# Run:  julia --project=. test/run_baseline.jl
#
# Establishes the euro/m2 a hand-set climate program earns on the twin, and the
# cost breakdown the RL agent will try to beat.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params = load_params(climate_json, crop_json)
params = update_params(params, ["psi2", "J_max"], [27800.0, 1.15e-4])   # nu4/k_ground now JSON defaults

gp = load_growth_params(growth_json)
gp = update_params(gp,
        ["node_rate", "veg_sink_max", "LAI_max", "set_start_day", "set_rate", "Wf_max"],
        [0.091,        21.0,           2.5,       16.0,            1.5,        16.0])

start_date = DateTime(1998, 7, 11, 0, 0)
ndays      = 100
weather    = load_weather(meteo, start_date, start_date + Day(ndays + 5); sky_from_clouds = true)

sp = daynight_setpoints(
    Tset_day = 22 + 273.15, Tset_night = 19 + 273.15,
    CO2_set  = 1200.0,
    light_start = 4.0, light_end = 20.0,
    VentpBand = 4.0, ofset = 1.0, ToutMax = 12 + 273.15,
)

@time res = simulate_twin(params, gp, weather, sp, datetime2unix(start_date), ndays; LAI0 = 0.5)

println("\n--- season summary ---")
println("final yield   = ", round(res.yield_FW[end], digits = 2), " kg FW/m2")
println("heating       = ", round(sum(res.heat_kWh), digits = 1), " kWh/m2")
println("CO2 dosed     = ", round(sum(res.co2_kg), digits = 1), " kg/m2")
println("lamp elec     = ", round(sum(res.lamp_kWh_peak) + sum(res.lamp_kWh_off), digits = 1), " kWh/m2")

e = season_economics(res; cropdays = ndays)
print_economics(e; label = "baseline fixed setpoints")
