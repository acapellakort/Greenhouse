# =============================================================================
# recalibrate_stomata.jl -- calibration-neutral fit of Jarvis gs_max.
#
# Keeps the literature response shapes (gs_D0, gs_I_half, gs_Ca_scale) and tunes
# ONLY gs_max so the responsive twin reproduces the FROZEN twin's seasonal gross
# assimilation at AMBIENT CO2 (400 ppm) -- where there is no enrichment closure,
# so this pins intrinsic canopy productivity. VPD/CO2 responsiveness then become
# correct additions on top rather than a recalibration of productivity.
#
#   Run:  julia -t auto --project=. test/recalibrate_stomata.jl
# =============================================================================
using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Random, Printf, Statistics

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

base = update_params(load_params(climate_json, crop_json), ["psi2","J_max"], [27800.0, 1.15e-4])
gp   = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max","leaf_per_node","leaf_lifespan_dd"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0, 0.040, 600.0])
bank = WeatherBank(meteo; plant_month = 8, plant_day = 14, ndays = 100, sky_from_clouds = true)

const YEARS = [1995, 2001, 2007, 2013]
const DENS  = 2.40
const WARM  = [16.01, 15.81, 782.73, 4.00, 4.26, 2.02, 11.75, 0.21, 52.11, 3.98, 0.81]

"Mean-over-years seasonal totals (gross assimilate g CH2O/m2, and net profit) for a program `a`."
function season_totals(p, a; years = YEARS)
    env = TwinEnv(p, gp, bank; forecast_skill = 0.0, rng = MersenneTwister(1), crop = :cohort)
    Pg = 0.0; prof = 0.0
    for y in years
        env_reset!(env; density_sched = (_ -> DENS), year = y)
        done = false; info = nothing
        while !done
            _, r, done, info = env_step!(env, a); prof += r; Pg += info.Pg
        end
    end
    return Pg / length(years), prof / length(years)
end

function main()
    p0 = base                                              # frozen (stomata_model = 0)
    warm_a = (WARM .- ACT_LO) ./ (ACT_HI .- ACT_LO)
    amb_a  = copy(warm_a); amb_a[3] = (400.0 - ACT_LO[3]) / (ACT_HI[3] - ACT_LO[3])  # CO2 -> 400 ppm

    target_Pg, _ = season_totals(p0, amb_a)                # frozen productivity @ ambient CO2
    @printf("frozen twin, ambient-CO2 reference: mean season Pg = %.1f g CH2O/m2\n", target_Pg)

    pj(gm) = update_params(base, ("stomata_model", "gs_max"), (1.0, gm))
    fPg(gm) = season_totals(pj(gm), amb_a)[1] - target_Pg  # increasing in gm

    lo, hi = 0.15, 1.5
    flo, fhi = fPg(lo), fPg(hi)
    @printf("bracket: Pg(gs_max=%.2f)-target = %+.1f ; Pg(gs_max=%.2f)-target = %+.1f\n", lo, flo, hi, fhi)
    if flo * fhi > 0
        println("!! no sign change in bracket -- widen [lo,hi] or check response shapes."); return
    end
    # bisection on gs_max (monotone increasing)
    gm = NaN
    for _ in 1:32
        gm = 0.5 * (lo + hi); fm = fPg(gm)
        fm > 0 ? (hi = gm) : (lo = gm)
        abs(fm) < 0.5 && break
    end
    @printf("\n>>> calibrated gs_max = %.4f  (was 0.33; matches frozen Pg within 0.5 g at ambient CO2)\n", gm)

    # report the isolated realism effect under the FULL CEM program (CO2=782)
    pgc, prc = season_totals(p0, warm_a)
    pgr, prr = season_totals(pj(gm), warm_a)
    println("\n--- CEM inc1 program (CO2 782 ppm), 4-year mean ---")
    @printf("  frozen         : Pg %.1f | profit %+.2f EUR/m2\n", pgc, prc)
    @printf("  responsive*    : Pg %.1f | profit %+.2f EUR/m2   (* recalibrated gs_max)\n", pgr, prr)
    @printf("  enrichment/VPD realism penalty on the OLD program: %+.2f EUR/m2\n", prr - prc)
    println("\nSet gs_max in constants_cropphoto.json to the calibrated value, then re-run CEM on the responsive twin.")
end

main()
