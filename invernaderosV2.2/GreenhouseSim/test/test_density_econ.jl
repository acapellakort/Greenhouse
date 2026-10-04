# =============================================================================
# test_density_econ.jl -- k_ext sensitivity + profit-proxy for the density choice
#   (layered cohort canopy on measured AiCU climate)
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

# economics coefficients (mirror V3.0 economics.jl)
PRICE=0.889; PLANT=0.22; CCO2=0.30
DENS = collect(2.0:0.5:5.0)

runq(dens, kext) = simulate_cohort_measured(params,
        update_params(gp,("stem_density","k_ext"),(dens,kext)), season; LAI0=0.5, layered=true)

# ---------------------------------------------- 1) k_ext x density (yield)
println("=== yield vs (k_ext, density): does higher extinction lower the optimum? ===")
@printf("%-6s |", "k_ext"); for d in DENS; @printf(" %5.1f", d); end; println("   | opt")
for k in (0.60, 0.75, 0.90, 1.10)
    ys = [runq(d, k).yield_FW[end] for d in DENS]
    @printf("%-6.2f |", k); for y in ys; @printf(" %5.1f", y); end
    @printf("   | %.1f stems\n", DENS[argmax(ys)])
end

# ---------------------------------------------- 2) PROFIT PROXY (base k=0.75)
# density-varying margin only: revenue - plant cost - CO2 dosing (~carbon fixed).
# heating/electricity/greenhouse-fixed are ~density-independent -> omitted (constant).
function proxy(dens, kext)
    r = runq(dens, kext)
    yld = r.yield_FW[end]
    co2fix = sum(r.Pg) * (44.0/30.0) / 1000.0          # kg CO2 / m2 fixed (dosing lower bound)
    rev = PRICE*yld; plants = PLANT*dens; co2c = CCO2*co2fix
    return (; dens, yld, co2fix, rev, plants, co2c, margin = rev - plants - co2c)
end
for (ktag, k) in (("k_ext=0.75 (base)",0.75), ("k_ext=0.90 (denser self-shades)",0.90))
    println("\n=== profit proxy, $ktag  [revenue - plants - CO2] ===")
    @printf("%-6s %6s %7s %8s %8s %8s\n","stems","yield","CO2kg","revenue","costs","margin")
    ms = Float64[]
    for d in DENS
        e = proxy(d, k); c = e.plants + e.co2c; push!(ms, e.margin)
        @printf("%-6.1f %6.1f %7.2f %8.2f %8.2f %8.2f\n", d, e.yld, e.co2fix, e.rev, c, e.margin)
    end
    @printf("-> profit-optimal density = %.1f stems (margin %.2f EUR/m2)\n", DENS[argmax(ms)], maximum(ms))
end
