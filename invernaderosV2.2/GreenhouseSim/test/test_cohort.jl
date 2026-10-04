# =============================================================================
# test_cohort.jl -- calibrate leaf production, then test de-leafing (cohort model)
# =============================================================================
using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using Printf

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA         = joinpath(INV, "invernaderos2.0", "data"); TEAM = "AiCU"
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")

params = update_params(load_params(climate_json, crop_json), ["J_max"], [1.15e-4])
gp0 = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0])
season = load_measured_climate(joinpath(DATA,TEAM,"Greenhouse_climate.csv"), joinpath(DATA,"meteo.csv"), params)

rold = simulate_growth_measured(params, gp0, season; LAI0=0.5)
@printf("OLD big-leaf reference:  yield %.1f kg  LAImax %.1f\n\n", rold.yield_FW[end], maximum(rold.LAI))

println("=== sweep leaf production (constant density 2.5, no de-leaf) ===")
@printf("%-8s %-8s  %6s %7s\n", "lpn", "lifesp", "LAImax", "yield")
for lpn in (0.03, 0.04, 0.05, 0.06), life in (400.0, 600.0)
    gp = update_params(gp0, ("leaf_per_node","leaf_lifespan_dd"), (lpn, life))
    r = simulate_cohort_measured(params, gp, season; LAI0=0.5, deleaf_target=Inf)
    @printf("%-8.3f %-8.0f  %6.1f %7.1f\n", lpn, life, maximum(r.LAI), r.yield_FW[end])
end

# pick a provisional calibration and test de-leafing on the realistic canopy
LPN, LIFE = 0.045, 500.0
gpc = update_params(gp0, ("leaf_per_node","leaf_lifespan_dd"), (LPN, LIFE))
@printf("\n=== de-leaf test at lpn=%.3f life=%.0f ===\n", LPN, LIFE)
function runc(tag, dens, deleaf)
    r = simulate_cohort_measured(params, update_params(gpc,("stem_density",),(dens,)), season; LAI0=0.5, deleaf_target=deleaf)
    @printf("%-32s yield %5.1f kg  LAImax %.1f\n", tag, r.yield_FW[end], maximum(r.LAI))
end
runc("constant 2.5",              2.5, Inf)
runc("constant 4.0",              4.0, Inf)
runc("dense 4.0 + deleaf->2.5",   4.0, 2.5)
runc("dense 4.0 + deleaf->3.0",   4.0, 3.0)
