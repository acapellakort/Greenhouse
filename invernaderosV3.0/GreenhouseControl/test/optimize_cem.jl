# =============================================================================
# optimize_cem.jl -- gradient-free (Cross-Entropy Method) baseline controller
#   Finds the best CONSTANT-season program (9 setpoints + planting density) to
#   maximize mean net profit across weather years. This is the "number to beat"
#   before SAC, and it audits the reward. Run threaded:  julia -t auto ...
# =============================================================================
using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates, Random, Printf, Statistics

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params = update_params(load_params(climate_json, crop_json), ["psi2","J_max"], [27800.0, 1.15e-4])
gp = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max","leaf_per_node","leaf_lifespan_dd"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0, 0.040, 600.0])
const CROP = :cohort   # optimize on the realistic (lamp-aware, layered) crop
bank = WeatherBank(meteo; plant_month = 8, plant_day = 14, ndays = 100, sky_from_clouds = true)

# --- CEM knobs ---------------------------------------------------------------
const DIM   = 12                       # 11 normalized setpoints + 1 normalized density
const POP   = 24
const NELITE= 6
const ITERS = 12
const YEARS = [1995, 2001, 2007, 2013] # fixed evaluation years (low-variance ranking)
const DENS_LO, DENS_HI = 1.5, 4.0

dens_of(x) = DENS_LO + clamp(x, 0, 1) * (DENS_HI - DENS_LO)

"Mean season net profit of a constant policy (first 9 dims) + density (10th) over YEARS."
function eval_policy(z; years = YEARS)
    a = z[1:11]; dens = dens_of(z[12])
    env = TwinEnv(params, gp, bank; forecast_skill = 0.0, rng = MersenneTwister(1), crop = CROP)
    tot = 0.0
    for y in years
        try
            env_reset!(env; density_sched = (_ -> dens), year = y)
            done = false
            while !done
                _, r, done, _ = env_step!(env, a)
                tot += r
            end
        catch
            return -1e3            # blew up -> heavy penalty
        end
    end
    return tot / length(years)
end

# --- run everything inside a function (avoids top-level soft-scope issues) ----
function main()
    base_phys = [22.0, 19.0, 1200.0, 16.0, 4.0, 1.0, 12.0, 0.5, 50.0, 4.0, 1.0]
    base_a    = (base_phys .- ACT_LO) ./ (ACT_HI .- ACT_LO)
    base_z    = vcat(base_a, (2.5 - DENS_LO)/(DENS_HI - DENS_LO))
    base_score = eval_policy(base_z)
    @printf("baseline program: mean profit over %s = %.2f EUR/m2\n\n", string(YEARS), base_score)

    rng = MersenneTwister(20260912)
    mu  = fill(0.5, DIM); sigma = fill(0.25, DIM)
    best_z = copy(mu); best_s = -Inf
    println("iter   mean-elite    best-so-far   (EUR/m2)")
    for it in 1:ITERS
        samples = [clamp.(mu .+ sigma .* randn(rng, DIM), 0.0, 1.0) for _ in 1:POP]
        scores  = fill(-Inf, POP)
        Threads.@threads for i in 1:POP
            scores[i] = eval_policy(samples[i])
        end
        order = sortperm(scores, rev = true)
        elite = samples[order[1:NELITE]]
        if scores[order[1]] > best_s
            best_s = scores[order[1]]; best_z = copy(samples[order[1]])
        end
        E = reduce(hcat, elite)
        mu    = vec(mean(E, dims = 2))
        sigma = vec(std(E, dims = 2)) .+ 1e-3
        @printf("%3d    %9.2f    %9.2f\n", it, mean(scores[order[1:NELITE]]), best_s)
    end

    bphys = ACT_LO .+ best_z[1:11] .* (ACT_HI .- ACT_LO)
    println("\n=== best CEM program ===")
    for (nm, v) in zip(ACT_NAMES, bphys); @printf("  %-12s %.2f\n", nm, v); end
    @printf("  %-12s %.2f stems/m2\n", "density", dens_of(best_z[12]))
    @printf("\nmean profit: baseline %.2f  ->  CEM %.2f  EUR/m2  (Δ %+.2f)\n",
            base_score, best_s, best_s - base_score)

    println("\nper-year profit (CEM best):")
    env = TwinEnv(params, gp, bank; forecast_skill = 0.0, crop = CROP)
    for y in YEARS
        env_reset!(env; density_sched = (_ -> dens_of(best_z[12])), year = y)
        tot = 0.0; done = false
        while !done; _, r, done, _ = env_step!(env, best_z[1:11]); tot += r; end
        @printf("  %d:  %.2f EUR/m2  (yield %.1f kg)\n", y, tot, env.gs.yield_FW)
    end
end

main()
