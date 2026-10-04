# =============================================================================
# sweep_leakage.jl  --  find the most physical leakage the corrected twin allows
# =============================================================================
# Run:  julia --project=. test/sweep_leakage.jl
#
# Context: nu4/nu4CO2 (heat & CO2 leakage) were inflated ~100x above Vanthoor
# (~1e-4) to 1e-2 for ODE stability, and the twin overrides them to 1e-3. That
# inflation was compensating for TWO bugs now fixed in V2.2 climate.jl:
#   (1) the ventilation-rate sign,   (2) the missing roof-vent term in h7.
# With both fixed, the greenhouse can shed heat through the roof and hold CO2, so
# the leakage should no longer need to be inflated. This sweep lowers nu4/nu4CO2
# from 1e-3 toward the physical 1e-4 and reports whether the twin stays stable,
# stops undershooting at night, and holds CO2 closer to setpoint.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Plots, Dates, Printf

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

# --- base parameters (same as run_twin, minus the leakage we are sweeping) ----
base_params = load_params(climate_json, crop_json)
base_params = update_params(base_params, ["psi2", "J_max"], [27800.0, 1.15e-4])

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

# --- leakage values to test (mg/current units) -------------------------------
leaks = [1.0e-3, 5.0e-4, 3.0e-4, 2.0e-4, 1.0e-4]

CO2_SET   = 1200.0
TSET_NIGHT = 19.0    # C, for night-undershoot reference

results = NamedTuple[]
pT = plot(title = "Daily mean air temp vs leakage", xlabel = "day", ylabel = "T air (C)")
pC = plot(title = "Daily mean CO2 vs leakage",       xlabel = "day", ylabel = "CO2 (mg/m3)")
hline!(pC, [CO2_SET], lc = :black, ls = :dash, label = "setpoint 1200")

println(@sprintf("%-9s %8s %8s %8s %9s %9s %8s %7s",
        "nu4", "Tmax", "Tmin", "CO2mean", "CO2min", "CO2max", "yield", "stable"))
for lk in leaks
    p  = update_params(base_params, ["nu4", "nu4CO2"], [lk, lk])
    ok = true
    local res
    try
        res = simulate_twin(p, gp, weather, sp, datetime2unix(start_date), ndays; LAI0 = 0.5)
        ok  = all(isfinite, res.Tmean) && all(isfinite, res.CO2mean) &&
              maximum(res.Tmean) < 60.0 && maximum(res.CO2mean) < 5000.0
    catch e
        ok = false
        @warn "twin failed at nu4=$lk" exception = e
    end
    if ok
        Tmax = maximum(res.Tmean); Tmin = minimum(res.Tmean)
        Cme  = sum(res.CO2mean)/length(res.CO2mean)
        Cmin = minimum(res.CO2mean); Cmax = maximum(res.CO2mean)
        yld  = res.yield_FW[end]
        push!(results, (; leak = lk, Tmax, Tmin, Cmean = Cme, Cmin, Cmax, yield = yld, stable = true))
        println(@sprintf("%-9.1e %8.1f %8.1f %8.0f %9.0f %9.0f %8.2f %7s",
                lk, Tmax, Tmin, Cme, Cmin, Cmax, yld, "yes"))
        lbl = @sprintf("nu4=%.0e", lk)
        plot!(pT, res.day, res.Tmean, lw = 2, label = lbl)
        plot!(pC, res.day, res.CO2mean, lw = 2, label = lbl)
    else
        push!(results, (; leak = lk, Tmax = NaN, Tmin = NaN, Cmean = NaN,
                          Cmin = NaN, Cmax = NaN, yield = NaN, stable = false))
        println(@sprintf("%-9.1e %8s %8s %8s %9s %9s %8s %7s",
                lk, "-", "-", "-", "-", "-", "-", "NO"))
    end
end

# reference lines on the temp plot
hline!(pT, [22.0], lc = :gray, ls = :dash, label = "day setpoint 22")
hline!(pT, [TSET_NIGHT], lc = :gray, ls = :dot, label = "night setpoint 19")

figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
savefig(plot(pT, pC, layout = (2, 1), size = (1000, 750)),
        joinpath(figdir, "sweep_leakage.png"))
println("\nFigure: ", joinpath(figdir, "sweep_leakage.png"))

println("""

Reading the result:
  * stable = yes AND CO2mean near 1200 AND Tmin not far below 19  => that leakage is fine.
  * As nu4 drops we expect CO2mean to RISE toward 1200 (less leak-out of dosed CO2)
    and the night Tmin to RISE toward 19 (less heat lost) -- both improvements.
  * Pick the SMALLEST (most physical) nu4 that is still 'stable = yes' and whose
    daytime Tmax is not driven up by lost leakage ventilation. That becomes the twin
    default and the starting point / prior center for the venting-day re-inference.
""")
