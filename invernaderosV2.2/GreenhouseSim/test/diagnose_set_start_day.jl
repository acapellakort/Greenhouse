# =============================================================================
# diagnose_set_start_day.jl  --  sensitivity of early yield to set_start_day
# =============================================================================
# Loads the posterior chain, fixes ALL parameters at their post-burn-in means
# EXCEPT set_start_day, and sweeps it over {12, 13, 14, 15, 16} days.
# Plots cumulative yield and weekly incremental yield so we can see which value
# best matches the early-season ramp (first 3-4 weekly bins).
#
# Run:  julia --project=invernaderosV2.2/GreenhouseSim \
#              invernaderosV2.2/GreenhouseSim/test/diagnose_set_start_day.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, JLD2, Plots, Statistics, Printf

# --- paths -------------------------------------------------------------------
INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA         = joinpath(INV, "invernaderos2.0", "data")
TEAM         = "AiCU"
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")
gh_csv       = joinpath(DATA, TEAM, "Greenhouse_climate.csv")
meteo_csv    = joinpath(DATA, "meteo.csv")
prod_csv     = joinpath(DATA, TEAM, "Production.csv")
cm_csv       = joinpath(DATA, TEAM, "CropManagement.csv")
chain_path   = joinpath(@__DIR__, "posteriors", TEAM * "_growth_chain.jld2")
figdir       = joinpath(@__DIR__, "figures"); mkpath(figdir)

# --- load --------------------------------------------------------------------
params  = load_params(climate_json, crop_json)
gp_base = load_growth_params(growth_json)

println("Loading measured climate ($TEAM) ...")
season = load_measured_climate(gh_csv, meteo_csv, params)

days_obs, y_obs = load_production(prod_csv)
cm_days, cm_Nleaves, cm_LAI = load_cropmanagement(cm_csv)

cum_at(days, cum, d) = (d <= days[1] ? float(cum[1]) :
                        d >= days[end] ? float(cum[end]) :
                        let i = searchsortedlast(days, d)
                            days[i] == d ? float(cum[i]) :
                            cum[i] + (d - days[i]) / (days[i+1] - days[i]) * (cum[i+1] - cum[i])
                        end)
weekly_incr(days, cum, edges) = diff([cum_at(days, cum, e) for e in edges])
edges  = collect(0:7:days_obs[end])
Δy_obs = weekly_incr(days_obs, y_obs, edges)

# --- load posterior mean (all params except set_start_day) -------------------
println("Loading chain: ", chain_path)
JLD2.@load chain_path chain QoI burn_in chain_length

QoI_growth = ["node_rate", "veg_sink_max", "LAI_max", "set_rate", "set_start_day"]
QoI_sigma  = ["sigma_y", "sigma_leaf", "sigma_lai"]
ng = length(QoI_growth)
@assert QoI[1:ng] == QoI_growth "Chain QoI order mismatch"

post      = chain[burn_in:end, :]
post_mean = vec(mean(post, dims = 1))

println("\nPosterior means (base for sweep):")
for (i, name) in enumerate(QoI)
    @printf "  %-16s = %.4f\n" name post_mean[i]
end

idx_ssd = findfirst(==("set_start_day"), QoI)

# --- sweep -------------------------------------------------------------------
sweep_days = [12.0, 13.0, 14.0, 15.0, 16.0]
colors     = [:steelblue, :darkorange, :forestgreen, :crimson, :purple]

# Residual sum of squares on weekly yield for each candidate
rss = Float64[]

println("\nSweep results (weekly yield RSS, full season):")
println("  set_start_day │  yield RSS  │  early RSS (wk 1-4)")
println("  ──────────────┼─────────────┼────────────────────")

pC = plot(days_obs, y_obs; lw = 2, lc = :black, label = "observed",
          xlabel = "day of season", ylabel = "kg FW / m²",
          title  = "Cumulative yield — set_start_day sweep ($TEAM)",
          legend = :topleft)

pW = plot(edges[2:end], Δy_obs; seriestype = :bar, fillalpha = 0.35,
          lc = :black, fc = :gray, label = "observed",
          xlabel = "day of season", ylabel = "kg FW / m² / wk",
          title  = "Weekly yield — set_start_day sweep ($TEAM)",
          legend = :topleft)

# zoom window for early season
early_end = min(35, Int(edges[end]))   # first 5 weeks

pE = plot(days_obs[days_obs .<= early_end], y_obs[days_obs .<= early_end];
          lw = 2, lc = :black, label = "observed",
          xlabel = "day of season", ylabel = "kg FW / m²",
          title  = "Early-season cumulative yield (days 0–$early_end)",
          legend = :topleft)

for (sd, col) in zip(sweep_days, colors)
    z_sweep      = copy(post_mean)
    z_sweep[idx_ssd] = sd
    gp_sweep = update_params(gp_base, QoI_growth, z_sweep[1:ng])
    res      = simulate_growth_measured(params, gp_sweep, season; LAI0 = 0.5)

    Δy_mod  = weekly_incr(res.day, res.yield_FW, edges)
    ry      = Δy_obs .- Δy_mod
    rss_all = sum(ry .^ 2)
    rss_early = sum(ry[1:min(4, length(ry))] .^ 2)
    push!(rss, rss_all)

    @printf "  %13.1f │  %9.4f  │  %.4f\n" sd rss_all rss_early

    lbl = "ssd = $(Int(sd))"
    plot!(pC, res.day, res.yield_FW;    lw = 2, lc = col, label = lbl)
    plot!(pW, edges[2:end], Δy_mod;     lw = 2, lc = col, label = lbl)
    plot!(pE, res.day[res.day .<= early_end],
              res.yield_FW[res.day .<= early_end]; lw = 2, lc = col, label = lbl)
end

savefig(pC, joinpath(figdir, "diag_ssd_cumulative.png"))
savefig(pW, joinpath(figdir, "diag_ssd_weekly.png"))
savefig(pE, joinpath(figdir, "diag_ssd_early.png"))
println("\nSaved: diag_ssd_cumulative.png, diag_ssd_weekly.png, diag_ssd_early.png")

# --- best value --------------------------------------------------------------
best_idx = argmin(rss)
best_ssd = sweep_days[best_idx]
println("\n=== Best set_start_day by weekly yield RSS: ", Int(best_ssd), " days ===")
println("    RSS values: ", [round(r, digits=4) for r in rss])
println("""
To fix set_start_day = $(Int(best_ssd)) in calibrate_growth.jl:
  1. Remove "set_start_day" from QoI_growth
  2. Add:  gp_base = update_params(gp_base, ["set_start_day"], [$(Int(best_ssd)).0])
  3. Remove the r_set likelihood term (no longer needed)
  4. n drops from 8 → 7 parameters
""")
