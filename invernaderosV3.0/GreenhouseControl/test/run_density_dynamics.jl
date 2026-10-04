# =============================================================================
# run_density_dynamics.jl -- LAI, photosynthesis, yield over time vs density
# =============================================================================
# Run:  julia --project=. test/run_density_dynamics.jl
# Shows the canopy-closure and light-capture dynamics behind the density effect.

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates, Plots

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

# canopy interception fraction from LAI (Beer's law, k = gp.k_ext)
intercept(lai) = 1 - exp(-gp.k_ext * lai)

scen = [("2.0", density_schedule([(0,2.0)])),
        ("2.5", density_schedule([(0,2.5)])),
        ("3.5", density_schedule([(0,3.5)])),
        ("5.0", density_schedule([(0,5.0)]))]

pL = plot(title="LAI", xlabel="day", ylabel="LAI"); pI = plot(title="light interception %", xlabel="day")
pY = plot(title="cumulative yield", xlabel="day", ylabel="kg/m2"); pP = plot(title="daily Pg", xlabel="day", ylabel="g CH2O/m2/d")
for (tag, ds) in scen
    res = simulate_twin(params, gp, weather, sp, datetime2unix(start_date), ndays; LAI0=0.5, density_sched=ds)
    plot!(pL, res.day, res.LAI, lw=2, label="dens $tag")
    plot!(pI, res.day, 100 .* intercept.(res.LAI), lw=2, label="dens $tag")
    plot!(pY, res.day, res.yield_FW, lw=2, label="dens $tag (end $(round(res.yield_FW[end],digits=1)))")
    plot!(pP, res.day, res.Pg, lw=2, label="dens $tag")
end
figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
savefig(plot(pL, pI, pY, pP, layout=(2,2), size=(1050,780)), joinpath(figdir, "density_dynamics.png"))
println("Figure: ", joinpath(figdir, "density_dynamics.png"))
