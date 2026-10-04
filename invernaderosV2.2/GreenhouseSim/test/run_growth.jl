# =============================================================================
# run_growth.jl  --  demo of the coupled climate + FvCB + cucumber growth twin
# =============================================================================
# Run from the GreenhouseSim folder:
#     julia --project=. test/run_growth.jl
#
# Runs a ~100-day cucumber season: each day integrates the fast climate+FvCB ODE
# at the current LAI to get gross assimilate, then advances the Marcelis-style
# source-sink growth model. Writes LAI / biomass / yield / assimilate plots.
#
# NOTE (v0): weather is the 1998 meteo series and controls are the 1-day
# observed set held constant (extrapolated) -- fine to demonstrate the coupled
# dynamics. Parameters in configfiles/cucumber_growth.json are literature
# placeholders; calibrate set_rate / veg_sink_max / Wf_max / SLA / asrq against
# the AGC-2018 Production.csv before trusting absolute yields.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using Dates, Plots

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")
obs_path     = joinpath(INV, "invernaderos2.0", "observed_data.csv")
growth_json  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")

# --- inputs ------------------------------------------------------------------
start_date = DateTime(1998, 7, 11, 0, 0)
ndays      = 100
end_date   = start_date + Day(ndays + 5)

params    = load_params(climate_json, crop_json)
gp        = load_growth_params(growth_json)
weather   = load_weather(meteo, start_date, end_date)
controls  = load_measured_controls(obs_path)

# --- run the coupled twin ----------------------------------------------------
@time res = simulate_growth(params, gp, weather, controls,
                            datetime2unix(start_date), ndays; LAI0 = 0.5)

# --- report ------------------------------------------------------------------
println("\n--- season summary (", ndays, " days) ---")
println("final LAI            = ", round(res.LAI[end], digits = 2))
println("final leaf DW        = ", round(res.W_leaf[end], digits = 1), " g/m2")
println("final stem DW        = ", round(res.W_stem[end], digits = 1), " g/m2")
println("cumulative yield     = ", round(res.yield_FW[end], digits = 2), " kg FW/m2")
println("mean daily Pg        = ", round(sum(res.Pg)/length(res.Pg), digits = 2), " g CH2O/m2/day")
println("AGC-2018 reference   ~ 36 kg FW/m2 over the full season (calibration target)")

# --- plots -------------------------------------------------------------------
figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)

pL = plot(res.day, res.LAI, xlabel = "day", ylabel = "LAI", label = "",
          title = "Canopy LAI")
pB = plot(res.day, res.W_leaf, xlabel = "day", ylabel = "g DW / m2",
          label = "leaf", title = "Biomass pools")
plot!(pB, res.day, res.W_stem,  label = "stem")
plot!(pB, res.day, res.W_fruit, label = "fruit on plant")
pY = plot(res.day, res.yield_FW, xlabel = "day", ylabel = "kg FW / m2",
          label = "model", title = "Cumulative fruit yield")
pP = plot(res.day, res.Pg, xlabel = "day", ylabel = "g CH2O / m2 / day",
          label = "", title = "Daily gross assimilate")

savefig(plot(pL, pB, pY, pP, layout = (2, 2), size = (1000, 700)),
        joinpath(figdir, "growth_season.png"))
savefig(pY, joinpath(figdir, "growth_yield.png"))
println("\nFigures written to ", figdir)
