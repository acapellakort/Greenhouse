# =============================================================================
# test_density.jl -- LAYERED canopy: big-leaf vs layered, LAI optimum,
#                     de-leaf (leaf only) and density schedule ("separate plants")
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
gp = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","set_start_day","set_rate","Wf_max","leaf_per_node","leaf_lifespan_dd"],
        [0.091, 21.0, 16.0, 1.5, 16.0, 0.040, 600.0])
season = load_measured_climate(joinpath(DATA,TEAM,"Greenhouse_climate.csv"), joinpath(DATA,"meteo.csv"), params)

runq(dens; sched=nothing, deleaf=Inf, layered=true) = simulate_cohort_measured(params,
        update_params(gp,("stem_density",),(dens,)), season; LAI0=0.5,
        density_sched=sched, deleaf_target=deleaf, layered=layered)

# ---------------------------------------------- 0) big-leaf vs layered Pg
println("=== big-leaf vs layered total Pg (same canopy) ===")
@printf("%-8s  %-9s %6s %8s\n", "stems", "Pg-model", "LAImax", "yield")
for dens in (2.5, 4.0), lay in (false, true)
    r = runq(dens; layered=lay)
    @printf("%-8.1f  %-9s %6.1f %8.1f\n", dens, lay ? "layered" : "big-leaf", maximum(r.LAI), r.yield_FW[end])
end

# ---------------------------------------------- 1) LAI optimum (layered)
println("\n=== constant-density sweep (LAYERED): where does denser stop paying? ===")
@printf("%-8s  %6s %7s\n", "stems", "LAImax", "yield")
rows = [(dens, runq(dens)) for dens in 2.0:0.5:5.0]
for (dens, r) in rows
    @printf("%-8.1f  %6.1f %7.1f\n", dens, maximum(r.LAI), r.yield_FW[end])
end
by = argmax([r.yield_FW[end] for (_, r) in rows])
@printf("-> best constant density = %.1f stems (%.1f kg)\n", rows[by][1], rows[by][2].yield_FW[end])

# ---------------------------------------------- 2) de-leaf (LEAF ONLY, no fruit loss)
println("\n=== de-leaf on the layered canopy (stems fixed, only shaded leaf removed) ===")
@printf("%-24s %6s %7s\n", "case", "LAImax", "yield")
for (tag, dens, dl) in (("constant 2.5",2.5,Inf),("constant 4.0",4.0,Inf),
                        ("4.0 deleaf->3.0",4.0,3.0),("4.0 deleaf->2.5",4.0,2.5),("4.0 deleaf->2.0",4.0,2.0))
    r = runq(dens; deleaf=dl)
    @printf("%-24s %6.1f %7.1f\n", tag, maximum(r.LAI), r.yield_FW[end])
end

# ---------------------------------------------- 3) density schedule (separate plants)
println("\n=== density schedule (LAYERED): start 4.0, separate to d_lo on day sep ===")
DHI = 4.0
r25 = runq(2.5); y25 = r25.yield_FW[end]
@printf("%-6s %-6s  %6s %7s   %s\n", "sep", "d_lo", "LAImax", "yield", "vs const-2.5")
@printf("%-6s %-6.1f  %6.1f %7.1f   (reference)\n", "const", 2.5, maximum(r25.LAI), y25)
rhi=runq(DHI); @printf("%-6s %-6.1f  %6.1f %7.1f   %+.1f\n","const",DHI,maximum(rhi.LAI),rhi.yield_FW[end],rhi.yield_FW[end]-y25)
for sep in (35,45,55), dlo in (2.0,2.5,3.0)
    r = runq(DHI; sched = d -> d < sep ? DHI : dlo)
    @printf("%-6d %-6.1f  %6.1f %7.1f   %+.1f\n", sep, dlo, maximum(r.LAI), r.yield_FW[end], r.yield_FW[end]-y25)
end
