# =============================================================================
# run_twin.jl  --  Stage-2 demo: full closed-loop cucumber digital twin
# =============================================================================
# Run:  julia --project=. test/run_twin.jl
#
# Runs a ~100-day cucumber season where the PID climate computer controls the
# greenhouse to a setpoint program AND the cucumber growth model rides on top:
# each day's controlled climate drives photosynthesis, the crop grows, and its
# LAI feeds back into the next day's climate. This is the environment the RL
# agent will act on (setpoints in, yield + resource use out).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Plots, Dates

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

# --- parameters --------------------------------------------------------------
params = load_params(climate_json, crop_json)
# actuator/physics + calibrated-source overrides for the control twin (V2.2 frozen):
params = update_params(params,
            ["psi2", "nu4", "nu4CO2", "J_max"],
            [27800.0, 1.0e-4, 1.0e-4, 1.15e-4])       # CO2 capacity, PHYSICAL Vanthoor leak (sweep_leakage.jl), calibrated Jmax

gp = load_growth_params(growth_json)
gp = update_params(gp,
        ["node_rate", "veg_sink_max", "LAI_max", "set_start_day", "set_rate", "Wf_max"],
        [0.091,        21.0,           2.5,       16.0,            1.5,        16.0])   # AiCU calibration

# --- season & weather --------------------------------------------------------
start_date = DateTime(1998, 7, 11, 0, 0)
ndays      = 100
weather    = load_weather(meteo, start_date, start_date + Day(ndays + 5); sky_from_clouds = true)

# --- setpoint program (cucumber: warm, CO2-enriched) -------------------------
sp = daynight_setpoints(
    Tset_day = 22 + 273.15, Tset_night = 19 + 273.15,
    CO2_set  = 1200.0,
    light_start = 4.0, light_end = 20.0,
    VentpBand = 4.0, ofset = 1.0, ToutMax = 12 + 273.15,
)

# --- run the coupled twin ----------------------------------------------------
@time res = simulate_twin(params, gp, weather, sp, datetime2unix(start_date), ndays; LAI0 = 0.5)

# --- report ------------------------------------------------------------------
println("\n--- twin season summary (", ndays, " days) ---")
println("final LAI        = ", round(res.LAI[end], digits = 2))
println("final yield      = ", round(res.yield_FW[end], digits = 2), " kg FW/m2")
println("final leaves     = ", round(res.node[end]))
println("mean air temp    = ", round(sum(res.Tmean)/length(res.Tmean), digits = 1), " C")
println("mean CO2         = ", round(sum(res.CO2mean)/length(res.CO2mean), digits = 0), " mg/m3")
println("mean daily Pg    = ", round(sum(res.Pg)/length(res.Pg), digits = 2), " g CH2O/m2/day")

# --- plots -------------------------------------------------------------------
figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
pL = plot(res.day, res.LAI, lw = 2, label = "", xlabel = "day", ylabel = "LAI", title = "Canopy LAI")
pY = plot(res.day, res.yield_FW, lw = 2, label = "", xlabel = "day", ylabel = "kg FW/m2",
          title = "Cumulative yield")
pB = plot(res.day, res.W_leaf, lw = 2, label = "leaf", xlabel = "day", ylabel = "g DW/m2",
          title = "Biomass")
plot!(pB, res.day, res.W_fruit, lw = 2, label = "fruit on plant")
pC = plot(res.day, res.Tmean, lw = 2, label = "T air (C)", xlabel = "day", title = "Daily climate")
plot!(pC, res.day, res.CO2mean ./ 100, lw = 2, label = "CO2 (x100 mg/m3)")
plot!(pC, res.day, res.Pg, lw = 2, label = "Pg (g CH2O/m2/d)")
savefig(plot(pL, pY, pB, pC, layout = (2, 2), size = (1000, 750)),
        joinpath(figdir, "twin_season.png"))

println("Figure: ", joinpath(figdir, "twin_season.png"))
