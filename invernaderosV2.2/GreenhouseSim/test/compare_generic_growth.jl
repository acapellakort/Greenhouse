# =============================================================================
# compare_generic_growth.jl -- legacy crop model vs Aaron's generic (Marcelis 1994) model (Note 66; forward)
# =============================================================================
# Run from the GreenhouseSim folder (~1 min, no MCMC):
#     julia --project=. test/compare_generic_growth.jl
#
# Same driver for every run: MEASURED AiCU air T, CO2 and canopy light (sun + lamps), days 0-114, as in the Note 21
# crop calibration. Shared crop parameters (node_rate, LAI_max, respiration, SLA, set gate, set_start_day 15) are
# the Note 21 posterior mean / JSON values, so the runs differ only in the partitioning model:
#   L    legacy V2.2: smoothstep fruit sink, 120 C d fruit period, capacity-fill vegetative sink (Note 21 fit)
#   G    generic, as ported (Aaron / M94 literature values): bell sink, 275 C d, harvest at stage 1,
#        sink and appearance = 0 when the daily mean T is outside 17.5-30 C
#   GW   G + harvest by size: a cohort is picked when its mean fruit reaches Wf_max = 16 g DM (= 400 g FW)
#   GWc  GW + T clamped to the M94 domain for sink and appearance (as M94's own equations)
#   In G, GW and GWc the node rate = flower-appearance rate (one flower per node); L keeps node_rate*(T - 10).
# Outputs: yield, weekly-yield fit, LAI fit, fruits set / harvested, harvested fruit size and age, and a figure.
# =============================================================================

ENV["OPTB_MODE"] = "setup"; ENV["N_DAYS"] = "115"
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

c21 = JLD2.load(N21_CHAIN); m21 = vec(mean(c21["chain"][c21["burn_in"]:end, 1:end-1], dims = 1))
g21 = [m21[findfirst(==(q), c21["QoI"])] for q in QoI_G]
gpL = update_params(gp_base, QoI_G, g21)
println("\nShared / legacy crop parameters (Note 21 mean): ", join(["$q=$(round(v, sigdigits=3))" for (q, v) in zip(QoI_G, g21)], "  "),
        "  stem_density=", gpL.stem_density, "  set_start_day=", gpL.set_start_day)
gpW  = update_params(gpL, ["harvest_W"], [gpL.Wf_max])
gpWc = update_params(gpW, ["gen_clampT"], [1.0])

function runcrop(gp, gen)
    s = init_growth_state(gp; LAI0 = 0.5)
    live = Dict{UInt,Tuple{Int,Fruit}}()
    day = Int[]; LAI = Float64[]; yld = Float64[]; node = Float64[]; res = Float64[]
    hn = Float64[]; hW = Float64[]; hage = Float64[]; nset = 0.0
    for d in 0:(N_DAYS - 1)
        M  = season_A.day_samples[d + 1]
        Pg = daily_assimilation_measured(M, s.LAI, params, season_A.dt, gp.k_ext)
        grow!(s, Pg, GreenhouseSim.daily_mean_T(M), gp, d; generic_params = gen)
        cur = Set(objectid(f) for f in s.fruits)
        for (k, (sd, f)) in collect(live)
            if !(k in cur); push!(hn, f.n); push!(hW, f.W); push!(hage, d - sd); delete!(live, k); end
        end
        for f in s.fruits
            k = objectid(f)
            if !haskey(live, k) && f.n > 0.0; live[k] = (d, f); nset += f.n; end
        end
        push!(day, d); push!(LAI, s.LAI); push!(yld, s.yield_FW); push!(node, s.node); push!(res, s.reserve)
    end
    nh = sum(hn; init = 0.0)
    return (; day, LAI, yield_FW = yld, node, reserve = res, nset, nharv = nh,
            fruitDW = nh > 0 ? sum(hW) / nh : NaN, age = nh > 0 ? sum(hn .* hage) / nh : NaN,
            onplant = isempty(s.fruits) ? 0.0 : sum(f.n for f in s.fruits))
end

runs = [("L    legacy (Note 21)", gpL, nothing),
        ("G    generic, harvest at stage 1", gpL, CUCUMBER_PARAMS),
        ("GW   generic, harvest at 16 g DM", gpW, CUCUMBER_PARAMS),
        ("GWc  GW + T clamped to M94 domain", gpWc, CUCUMBER_PARAMS)]
R = [runcrop(gp, gen) for (_, gp, gen) in runs]

Tday = [GreenhouseSim.daily_mean_T(season_A.day_samples[d + 1]) for d in 0:(N_DAYS - 1)]
@printf("\nDaily mean air T: mean %.1f C, min %.1f, max %.1f; days below 17.5 C (generic sink/appearance = 0): %d of %d\n",
        mean(Tday), minimum(Tday), maximum(Tday), count(<(17.5), Tday), N_DAYS)

lai_at(r, d) = r.LAI[min(d, N_DAYS - 1) + 1]
println("\n", "="^150)
@printf("  %-34s | %7s %8s | %8s %7s | %6s %6s | %8s %8s %7s | %8s %8s | %6s\n", "run", "yield", "d42-114", "wk RMSE", "wk bias",
        "LAI rm", "LAImean", "set /m2", "harv/m2", "onplant", "fruit DW", "fruit FW", "age d")
println("  ", "-"^148)
for ((nm, gp, _), r) in zip(runs, R)
    wk = weekly_incr(r.day, r.yield_FW, edges)
    larm = sqrt(mean(abs2, [lai_at(r, d) for d in cm_days] .- cm_L))
    @printf("  %-34s | %7.2f %8.2f | %8.3f %+7.3f | %6.2f %6.2f | %8.1f %8.1f %7.1f | %8.1f %8.0f | %6.1f\n", nm,
            r.yield_FW[end], r.yield_FW[end] - r.yield_FW[43], sqrt(mean(abs2, wk .- Δy_obs)), mean(wk .- Δy_obs),
            larm, mean(r.LAI), r.nset, r.nharv, r.onplant, r.fruitDW, r.fruitDW / gp.DMC_fruit, r.age)
end
yo = cum_at(days_obs, y_obs, N_DAYS - 1)
@printf("  %-34s | %7.2f %8.2f |\n", "observed (Production.csv)", yo, yo - cum_at(days_obs, y_obs, 42))
println("\n  Nodes vs observed leaf number (generic runs: node = flower appearance, Note 66)")
for ((nm, _, _), r) in zip(runs, R)
    nm_ = [r.node[d + 1] for d in cm_days]
    @printf("  %-34s | node d%d %6.1f (obs %5.1f) | RMSE %5.2f | bias %+5.2f\n", nm, cm_days[end], nm_[end], cm_N[end],
            sqrt(mean(abs2, nm_ .- cm_N)), mean(nm_ .- cm_N))
end

println("\n  Weekly yield (kg FW m-2 wk-1)")
@printf("  %8s | %7s | %7s %7s %7s %7s\n", "week end", "obs", "L", "G", "GW", "GWc")
W = [weekly_incr(r.day, r.yield_FW, edges) for r in R]
for i in eachindex(Δy_obs)
    @printf("  %8d | %7.2f | %7.2f %7.2f %7.2f %7.2f\n", edges[i + 1], Δy_obs[i], W[1][i], W[2][i], W[3][i], W[4][i])
end

figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
labs = ["L", "G", "GW", "GWc"]
pC = plot(days_obs, y_obs, lw = 3, color = :black, label = "observed", xlabel = "day", ylabel = "kg FW m-2", title = "Cumulative yield")
for (r, l) in zip(R, labs); plot!(pC, r.day, r.yield_FW, lw = 2, label = l); end
pW = scatter(edges[2:end], Δy_obs, color = :black, label = "observed", xlabel = "day", ylabel = "kg FW m-2 wk-1", title = "Weekly yield")
for (w, l) in zip(W, labs); plot!(pW, edges[2:end], w, lw = 2, label = l); end
pL = scatter(cm_days, cm_L, color = :black, label = "derived obs", xlabel = "day", ylabel = "LAI", title = "LAI")
for (r, l) in zip(R, labs); plot!(pL, r.day, r.LAI, lw = 2, label = l); end
pN = scatter(cm_days, cm_N, color = :black, label = "observed leaves", xlabel = "day", ylabel = "nodes per stem", title = "Nodes")
for (r, l) in zip(R, labs); plot!(pN, r.day, r.node, lw = 2, label = l); end
fig = plot(pC, pW, pL, pN, layout = (1, 4), size = (1900, 420), margin = 5Plots.mm)
fn = joinpath(figdir, "compare_generic_growth.png"); savefig(fig, fn); println("\n  figure: ", fn)
println("""
  Read: G is Aaron's model as ported. Its fruits keep growing for 275 C d toward ~60 g DM, so harvest by age gives
  very large fruits and few of them; GW picks at the Note 21 harvest standard (16 g DM). Literature values only --
  nothing here is fitted. If GW/GWc are in the right range, the next step is a generic-path calibration
  (long run; ask first).""")
