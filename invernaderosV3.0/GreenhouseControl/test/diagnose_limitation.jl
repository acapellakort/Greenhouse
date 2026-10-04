# =============================================================================
# diagnose_limitation.jl -- is the autumn scenario light-limited?
#
#   Part A (env level): season gross assimilate Pg & net profit under the CEM
#     program at CO2 setpoints 400/700/1000/1200 ppm, FROZEN model (flag 0).
#     Flat Pg across CO2 => not CO2-limited => light-limited.
#   Part B (leaf level): FvCB limitation split (Aj light / Ac RuBisCO / Ap TPU)
#     over an autumn daytime grid; which term is the min is the binding process.
#
#   Run:  julia -t auto --project=. test/diagnose_limitation.jl
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

function season_totals(p, a; years = YEARS)
    env = TwinEnv(p, gp, bank; forecast_skill = 0.0, rng = MersenneTwister(1), crop = :cohort)
    Pg = 0.0; prof = 0.0
    for y in years
        env_reset!(env; density_sched = (_ -> DENS), year = y); done = false; info = nothing
        while !done; _, r, done, info = env_step!(env, a); prof += r; Pg += info.Pg; end
    end
    return Pg / length(years), prof / length(years)
end

# ---- leaf-level FvCB limitation split (mirror of assimilation, frozen g_eff) ----
function fvcb_split(T, C, Iw, p)
    (; O_a, Sco25, E_Soc, Rgas, alpha, J_max, theta,
       V_cmax25, Q10_Vcmax, Rd_day, K_C25, Q10_KC, K_O25, Q10_KO) = p
    I = 4.6e-6*Iw; C_i = 0.509e-6*C
    Gamma_st = (0.5*O_a)/(Sco25*exp(((T-298.15)/T)*E_Soc/(298.15*Rgas)))
    J_t = ((alpha*I+J_max) - sqrt((alpha*I+J_max)^2 - 4*theta*alpha*I*J_max))/(2*theta)
    g = GreenhouseControl.GreenhouseSim.g_eff()
    V_cmax = V_cmax25*Q10_Vcmax^(0.1*(T-298.15))/(1+exp(0.128*(T-315.15)))
    K_C = K_C25*Q10_KC^(0.1*(T-298.15)); K_O = K_O25*Q10_KO^(0.1*(T-298.15))
    pc = -(V_cmax + g*(C_i+K_C*(1+O_a/K_O)) - Rd_day)
    qc = g*(V_cmax*max(C_i-Gamma_st,0) - (C_i+K_C*(1+O_a/K_O))*Rd_day)
    Ac = pc^2-4*qc >= 0 ? 0.5*(-pc-sqrt(pc^2-4*qc)) : -Rd_day
    pj = -0.25*J_t + Rd_day - g*(C_i+2*Gamma_st)
    qj = 0.25*g*max(C_i-Gamma_st,0)*J_t - g*(C_i+2*Gamma_st)*Rd_day
    Aj = pj^2-4*qj >= 0 ? 0.5*(-pj-sqrt(pj^2-4*qj)) : -Rd_day
    Ap = 1.5*V_cmax - Rd_day
    return Aj, Ac, Ap
end

function main()
    println("== Part A: season Pg & profit vs CO2 setpoint (FROZEN model) ==")
    @printf("%8s %12s %12s\n", "CO2 ppm", "Pg gCH2O/m2", "profit EUR/m2")
    for co2 in (400.0, 700.0, 1000.0, 1200.0)
        a = (WARM .- ACT_LO) ./ (ACT_HI .- ACT_LO)
        a[3] = (co2 - ACT_LO[3]) / (ACT_HI[3] - ACT_LO[3])
        pg, pr = season_totals(base, a)
        @printf("%8.0f %12.1f %12.2f\n", co2, pg, pr)
    end
    println("  (flat Pg across CO2 => not CO2-limited => light-limited)")

    println("\n== Part B: FvCB leaf limitation split (autumn daytime grid) ==")
    println("  min of {Aj light, Ac RuBisCO, Ap TPU} is the binding process; Iw = canopy-top light")
    @printf("  %6s %6s | %8s %8s %8s  limiting\n", "Iw", "CO2", "Aj", "Ac", "Ap")
    for Iw in (20.0, 60.0, 120.0, 250.0, 500.0), co2 in (400.0, 800.0)
        Aj, Ac, Ap = fvcb_split(293.0, co2*1.83, Iw, base)
        lim = argmin((Aj, Ac, Ap)); name = ("Aj-light","Ac-RuBisCO","Ap-TPU")[lim]
        @printf("  %6.0f %6.0f | %8.4f %8.4f %8.4f  %s\n", Iw, co2, Aj, Ac, Ap, name)
    end
    println("\n  Transmitted canopy-top light in autumn is mostly < ~150 W/m2 by day;")
    println("  if Aj is the min there, the crop is light-limited and stomata/CO2 are second-order.")
end

main()
