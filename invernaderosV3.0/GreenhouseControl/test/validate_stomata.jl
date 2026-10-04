# =============================================================================
# validate_stomata.jl -- checks the responsive-stomata / diurnal-respiration
# extension of the (frozen-by-default) V2.2 engine.
#
#   TEST 1  flag=0 is bit-for-bit: VPD is ignored and gs == g_eff.
#   TEST 2  Jarvis responds in the right direction (VPD down, light up, CO2 down).
#   TEST 3  anchor: gs at reference conditions ~= old constant g_eff (0.136).
#   TEST 4  end-to-end TwinEnv episode runs & stays finite with flag=0 and flag=1.
#
#   Run (after the campaign, to avoid CPU contention):
#     julia -t auto --project=. test/validate_stomata.jl
# =============================================================================
using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Random, Printf, Statistics

const GS = GreenhouseControl.GreenhouseSim   # photosynthesis internals live here

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

p0 = base                                                # stomata_model = 0 (frozen default)
p1 = update_params(base, ("stomata_model",), (1.0,))     # Jarvis on
mgm3(ppm) = ppm * 1.83                                   # ppm -> mg/m3 (~20C)

function main()
    # ---- TEST 1: flag=0 bit-for-bit --------------------------------------
    println("== TEST 1: constant-stomata (flag=0) is bit-for-bit ==")
    maxdiff = 0.0
    for T in (288.0, 298.0, 305.0), C in (600.0, 1200.0, 2000.0),
        Iw in (0.0, 50.0, 300.0), VPD in (0.0, 1.5, 3.0)
        a0 = GS.assimilation(T, C, Iw, 2.0, p0)
        av = GS.assimilation(T, C, Iw, 2.0, p0; VPD = VPD)
        maxdiff = max(maxdiff, abs(av - a0))
    end
    gc = GS.stomatal_conductance(150.0, mgm3(400), 0.8, 298.0, p0)
    ge = GS.g_eff()
    @printf("  max |A(VPD) - A(noVPD)| over grid = %.3e   (expect 0)\n", maxdiff)
    @printf("  gs(flag=0) = %.6f ; g_eff = %.6f ; identical = %s\n", gc, ge, gc == ge)
    println(maxdiff == 0.0 && gc == ge ? "  PASS" : "  FAIL")

    # ---- TEST 2: Jarvis monotonicity (flag=1) ----------------------------
    println("\n== TEST 2: Jarvis responds in the correct direction ==")
    jg(Iw, ppm, VPD) = GS.stomatal_conductance(Iw, mgm3(ppm), VPD, 298.0, p1)  # total (series)
    println("  VPD sweep (Iw=150, 400 ppm) -- expect DECREASING:")
    gv = [jg(150.0, 400.0, v) for v in (0.2, 0.8, 1.5, 3.0)]
    @printf("    VPD 0.2/0.8/1.5/3.0 kPa -> %.4f %.4f %.4f %.4f\n", gv...)
    println("  light sweep (VPD=0.8, 400 ppm) -- expect INCREASING:")
    gl = [jg(i, 400.0, 0.8) for i in (0.0, 30.0, 150.0, 600.0)]
    @printf("    Iw 0/30/150/600 W/m2 -> %.4f %.4f %.4f %.4f\n", gl...)
    println("  CO2 sweep (Iw=150, VPD=0.8) -- expect DECREASING:")
    gc2 = [jg(150.0, c, 0.8) for c in (400.0, 700.0, 1000.0, 1400.0)]
    @printf("    CO2 400/700/1000/1400 ppm -> %.4f %.4f %.4f %.4f\n", gc2...)
    ok2 = issorted(gv, rev=true) && issorted(gl) && issorted(gc2, rev=true)
    println(ok2 ? "  PASS" : "  FAIL")

    # ---- TEST 3: anchor to old calibration -------------------------------
    println("\n== TEST 3: reference-condition anchor ==")
    gref = jg(150.0, 400.0, 0.8)
    @printf("  gs(Iw=150, 400 ppm, VPD=0.8) = %.4f  vs  g_eff = %.4f  (ratio %.2f)\n",
            gref, ge, gref/ge)
    println(0.7 <= gref/ge <= 1.4 ? "  PASS (within ~30% of frozen operating point)" :
                                    "  WARN (anchor drifted; consider retuning gs_max)")

    # ---- TEST 4: end-to-end episode, flag 0 vs 1 -------------------------
    println("\n== TEST 4: TwinEnv episode (CEM inc1 program) ==")
    WARM = [16.01, 15.81, 782.73, 4.00, 4.26, 2.02, 11.75, 0.21, 52.11, 3.98, 0.81]
    a    = (WARM .- ACT_LO) ./ (ACT_HI .- ACT_LO)
    DENS = 2.40
    function run_ep(p, year)
        env = TwinEnv(p, gp, bank; forecast_skill = 0.0, rng = MersenneTwister(1), crop = :cohort)
        env_reset!(env; density_sched = (_ -> DENS), year = year)
        tot = 0.0; done = false; nfail = 0; info = nothing
        while !done
            _, r, done, info = env_step!(env, a); tot += r
            (haskey(info, :failed) && info.failed) && (nfail += 1)
        end
        return tot, env.gs.LAI, nfail
    end
    for y in (1995, 2001, 2007)
        t0, l0, f0 = run_ep(p0, y)
        t1, l1, f1 = run_ep(p1, y)
        @printf("  year %d | flag0 profit %7.2f LAI %.2f fails %d | flag1 profit %7.2f LAI %.2f fails %d\n",
                y, t0, l0, f0, t1, l1, f1)
    end
    println("  (flag0 should reproduce the frozen model; flag1 should run finite with 0 fails)")
end

main()
