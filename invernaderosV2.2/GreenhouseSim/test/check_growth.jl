# =============================================================================
# check_growth.jl  --  FAST visual check of the growth model (no MCMC)
# =============================================================================
# Run:  julia --project=. test/check_growth.jl
#
# Runs ONE season of the measured-climate growth model at the current best-fit
# point estimates and plots weekly/cumulative yield, leaf number and LAI against
# the AGC-2018 data. Use this to eyeball the new structure (harvest flushes from
# carbohydrate-regulated set, and the late-season topping) and to tune the flush
# parameters (set_r_low / set_r_high / set_rate) in configfiles/cucumber_growth.json
# quickly, BEFORE launching a full calibration.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Dates, Plots

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA         = joinpath(INV, "invernaderos2.0", "data")
TEAM         = "AiCU"
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")
gh_csv       = joinpath(DATA, TEAM, "Greenhouse_climate.csv")
meteo_csv    = joinpath(DATA, "meteo.csv")
prod_csv     = joinpath(DATA, TEAM, "Production.csv")
cm_csv       = joinpath(DATA, TEAM, "CropManagement.csv")

params = load_params(climate_json, crop_json)
gp     = load_growth_params(growth_json)

# current best-fit point estimates (from the calibration runs) for the visual;
# flush params (set_r_low/high, top_day) come from the JSON so you can tune them.
params = update_params(params, ["J_max"], [1.15e-4])
gp     = update_params(gp, ["node_rate", "veg_sink_max", "LAI_max", "set_start_day", "set_rate", "Wf_max"],
                            [0.091,        21.0,           2.5,       16.0,            1.5,        16.0])

season = load_measured_climate(gh_csv, meteo_csv, params)
res    = simulate_growth_measured(params, gp, season; LAI0 = 0.5)

# --- observed ----------------------------------------------------------------
days_obs, y_obs         = load_production(prod_csv)
cm_days, cm_Nleaves, cm_LAI = load_cropmanagement(cm_csv)

cum_at(days, cum, d) = (d <= days[1] ? float(cum[1]) :
                        d >= days[end] ? float(cum[end]) :
                        let i = searchsortedlast(days, d)
                            days[i] == d ? float(cum[i]) :
                            cum[i] + (d - days[i]) / (days[i+1] - days[i]) * (cum[i+1] - cum[i])
                        end)
weekly_incr(days, cum, edges) = diff([cum_at(days, cum, e) for e in edges])
edges = collect(0:7:days_obs[end])

figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)

pW = plot(edges[2:end], weekly_incr(days_obs, y_obs, edges), seriestype = :scatter, ms = 4,
          label = "observed", xlabel = "day", ylabel = "kg FW/m2/wk",
          title = "Weekly yield check ($TEAM)")
plot!(pW, edges[2:end], weekly_incr(res.day, res.yield_FW, edges), lw = 2, label = "model")
savefig(pW, joinpath(figdir, "check_weekly.png"))

pC = plot(days_obs, y_obs, lw = 2, label = "observed", xlabel = "day",
          ylabel = "kg FW/m2", title = "Cumulative yield check ($TEAM)")
plot!(pC, res.day, res.yield_FW, lw = 2, ls = :dash, label = "model")
savefig(pC, joinpath(figdir, "check_cumulative.png"))

pN = plot(cm_days, cm_Nleaves, seriestype = :scatter, ms = 4, label = "observed",
          xlabel = "day", ylabel = "leaves/stem", title = "Leaf number check ($TEAM)")
plot!(pN, res.day, res.node, lw = 2, label = "model")
savefig(pN, joinpath(figdir, "check_leaf.png"))

pL = plot(cm_days, cm_LAI, seriestype = :scatter, ms = 4, label = "observed (derived)",
          xlabel = "day", ylabel = "LAI", title = "LAI check ($TEAM)")
plot!(pL, res.day, res.LAI, lw = 2, label = "model")
savefig(pL, joinpath(figdir, "check_lai.png"))

println("final yield = ", round(res.yield_FW[end], digits = 2), " kg/m2 (obs ",
        round(y_obs[end], digits = 2), ");  final leaves = ", round(res.node[end]),
        " (obs ", round(cm_Nleaves[end]), ")")
println("Figures: check_weekly / check_cumulative / check_leaf / check_lai .png in ", figdir)
