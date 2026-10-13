# =============================================================================
# forecast_growth.jl  --  Probabilistic production forecast from posterior chain
# =============================================================================
# Loads the t-walk posterior for AiCU growth parameters, draws N_SAMPLES from
# the post-burn-in chain, and produces two kinds of output:
#
#   1. Full-season POSTERIOR PREDICTIVE (oracle climate) — model validation:
#      does the observed trajectory fall inside the 90% band?
#
#   2. Honest N-week FORECAST from FORECAST_DAY:
#      future greenhouse climate replaced by the within-season mean up to that
#      day (Level-1 / climatological approach). This avoids using measured
#      future T, CO2, and light that would be unknown at forecast time.
#      Both oracle and climatological bands are plotted so the cost of the
#      approximation is visible.
#
# NOTE on Level-1 assumption: future T_in, CO2_in, light ≈ mean of past
#   days 0:(FORECAST_DAY-1). This implicitly assumes the grower continues
#   their control program and outdoor weather is seasonally average.
#   A Level-2 forecast would propagate the climate ODE from FORECAST_DAY
#   using planned setpoints + outdoor weather forecast (future work).
#
# Run:  julia --project=invernaderosV2.2/GreenhouseSim \
#              invernaderosV2.2/GreenhouseSim/test/forecast_growth.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, JLD2, Plots, Statistics, Random, Printf

# =============================================================================
# CONFIG
# =============================================================================
N_SAMPLES    = 200     # posterior draws
FORECAST_DAY = 42      # day on which the forecast is "issued" (0-based season day)
FORECAST_WKS = 4       # weeks ahead

# =============================================================================
# paths
# =============================================================================
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

# =============================================================================
# load base params + data
# =============================================================================
params  = load_params(climate_json, crop_json)
gp_base = load_growth_params(growth_json)
gp_base = update_params(gp_base, ["set_start_day"], [15.0])  # must match calibrate_growth.jl

println("Loading measured climate (", TEAM, ") ...")
season = load_measured_climate(gh_csv, meteo_csv, params)
println("  season length = ", season.ndays, " days")

days_obs, y_obs          = load_production(prod_csv)
cm_days, cm_Nleaves, cm_LAI = load_cropmanagement(cm_csv)

# interpolation helpers (same as calibrate_growth.jl)
cum_at(days, cum, d) = (d <= days[1] ? float(cum[1]) :
                        d >= days[end] ? float(cum[end]) :
                        let i = searchsortedlast(days, d)
                            days[i] == d ? float(cum[i]) :
                            cum[i] + (d - days[i]) / (days[i+1] - days[i]) * (cum[i+1] - cum[i])
                        end)
at_day(v, d) = v[min(d, length(v) - 1) + 1]

# =============================================================================
# Level-1 climatological season
# For days >= forecast_day: replace day_samples with the mean climate of
# days 0:(forecast_day-1). Rows are [Tair_K; CO2_mg/m3; Iw].
# =============================================================================
function climatological_season(season::MeasuredSeason, forecast_day::Int)
    forecast_day >= season.ndays && return season   # nothing to replace

    pre_mats = filter(m -> size(m, 2) > 0, season.day_samples[1:forecast_day])
    isempty(pre_mats) && error("No valid samples before forecast_day=$forecast_day")

    all_pre = hcat(pre_mats...)                          # 3 × N_total
    μ       = vec(mean(all_pre, dims = 2))               # mean [T; CO2; Iw]

    # replicate to a typical day length (most common sample count in pre window)
    n_typ   = round(Int, mean(size(m, 2) for m in pre_mats))
    μ_day   = repeat(reshape(μ, 3, 1), 1, n_typ)        # 3 × n_typ constant day

    new_samples = copy(season.day_samples)
    for d in forecast_day:(season.ndays - 1)
        new_samples[d + 1] = μ_day
    end
    return MeasuredSeason(season.ndays, new_samples, season.dt)
end

season_clim = climatological_season(season, FORECAST_DAY)
println("Climatological season built: days 0–$(FORECAST_DAY-1) measured, ",
        "days $(FORECAST_DAY)–$(season.ndays-1) climatological mean.")

# =============================================================================
# load posterior chain
# =============================================================================
println("Loading chain: ", chain_path)
JLD2.@load chain_path chain QoI burn_in chain_length

QoI_growth = ["node_rate", "veg_sink_max", "LAI_max", "set_rate"]
# Wf_max        fixed at 16.0 g DW/fruit (400 g FW × DMC 0.04) — not inferred
# set_start_day fixed at 15 d (grower observation + RSS sweep)   — not inferred
ng = length(QoI_growth)
@assert QoI[1:ng] == QoI_growth "Chain QoI order mismatch — re-check calibrate_growth.jl"

post   = chain[burn_in:end, 1:ng]
n_post = size(post, 1)
println("  post-burn-in rows: ", n_post, " → drawing ", N_SAMPLES)

rng = MersenneTwister(42)
idx = N_SAMPLES <= n_post ?
      sort(randperm(rng, n_post)[1:N_SAMPLES]) :
      rand(rng, 1:n_post, N_SAMPLES)

# =============================================================================
# forward ensemble — oracle (full measured) AND climatological
# =============================================================================
day_axis = 0:season.ndays
nd       = length(day_axis)

yield_oracle = Matrix{Float64}(undef, nd, N_SAMPLES)
LAI_oracle   = Matrix{Float64}(undef, nd, N_SAMPLES)
yield_clim   = Matrix{Float64}(undef, nd, N_SAMPLES)
LAI_clim     = Matrix{Float64}(undef, nd, N_SAMPLES)

function run_ensemble!(yield_mat, LAI_mat, szn, post, idx, params, gp_base, QoI_growth, day_axis, label)
    println("Running $N_SAMPLES simulations ($label) ...")
    for (j, i) in enumerate(idx)
        gp  = update_params(gp_base, QoI_growth, post[i, :])
        res = simulate_growth_measured(params, gp, szn; LAI0 = 0.5)
        for (di, d) in enumerate(day_axis)
            yield_mat[di, j] = cum_at(res.day, res.yield_FW, d)
            LAI_mat[di, j]   = Float64(at_day(res.LAI, d))
        end
        j % 50 == 0 && println("  ", j, " / ", N_SAMPLES)
    end
    println("Done ($label).")
end

run_ensemble!(yield_oracle, LAI_oracle, season,      post, idx, params, gp_base, QoI_growth, day_axis, "oracle")
run_ensemble!(yield_clim,   LAI_clim,   season_clim, post, idx, params, gp_base, QoI_growth, day_axis, "climatological")

pct(mat, p) = [quantile(view(mat, d, :), p) for d in 1:nd]
days = collect(day_axis)

# oracle bands
# oracle bands
p05_yo = pct(yield_oracle, 0.05); p25_yo = pct(yield_oracle, 0.25)
p50_yo = pct(yield_oracle, 0.50)
p75_yo = pct(yield_oracle, 0.75); p95_yo = pct(yield_oracle, 0.95)
p05_Lo = pct(LAI_oracle,   0.05); p50_Lo = pct(LAI_oracle,   0.50); p95_Lo = pct(LAI_oracle,   0.95)

# climatological bands
p05_yc = pct(yield_clim, 0.05); p25_yc = pct(yield_clim, 0.25)
p50_yc = pct(yield_clim, 0.50)
p75_yc = pct(yield_clim, 0.75); p95_yc = pct(yield_clim, 0.95)
p05_Lc = pct(LAI_clim,   0.05); p50_Lc = pct(LAI_clim,   0.50); p95_Lc = pct(LAI_clim,   0.95)

# =============================================================================
# PLOT 1: full-season posterior predictive (oracle) — model validation
# =============================================================================
pY = plot(days, p50_yo; lw = 2, lc = :steelblue, label = "median (oracle)",
          xlabel = "day of season", ylabel = "kg FW / m²",
          title  = "Posterior predictive — cumulative yield ($TEAM)",
          legend = :topleft)
plot!(pY, days, p05_yo; fillrange = p95_yo,
      fillalpha = 0.20, lw = 0, fillcolor = :steelblue, label = "5–95 % (oracle)")
scatter!(pY, days_obs, y_obs; ms = 5, mc = :black, label = "observed")
vline!(pY, [FORECAST_DAY]; lw = 1.5, lc = :crimson, ls = :dash,
       label = "forecast day $FORECAST_DAY")
savefig(pY, joinpath(figdir, "forecast_yield_fullseason.png"))
println("Saved: forecast_yield_fullseason.png")

# =============================================================================
# PLOT 2: full-season LAI posterior predictive (oracle)
# =============================================================================
pL = plot(days, p50_Lo; lw = 2, lc = :forestgreen, label = "median (oracle)",
          xlabel = "day of season", ylabel = "LAI",
          title  = "Posterior predictive — LAI ($TEAM)", legend = :topleft)
plot!(pL, days, p05_Lo; fillrange = p95_Lo,
      fillalpha = 0.20, lw = 0, fillcolor = :forestgreen, label = "5–95 %")
scatter!(pL, cm_days, cm_LAI; ms = 5, mc = :black, label = "observed (derived)")
savefig(pL, joinpath(figdir, "forecast_LAI_fullseason.png"))
println("Saved: forecast_LAI_fullseason.png")

# =============================================================================
# PLOT 3: forecast window — oracle vs climatological bands
# =============================================================================
fd   = FORECAST_DAY
fend = min(fd + FORECAST_WKS * 7, season.ndays)
row(d) = d - first(day_axis) + 1
wd = collect(fd:fend); wr = row(fd):row(fend)
pF = plot(wd, p50_yo[wr]; lw = 2, lc = :steelblue, ls = :dash, label = "oracle median",
          xlabel = "day of season", ylabel = "kg FW / m²",
          title  = "$(FORECAST_WKS)-week forecast from day $fd ($TEAM)",
          legend = :topleft)
# oracle: 90 % band (light), 50 % band on top (same alpha → appears darker)
plot!(pF, wd, p05_yo[wr]; fillrange = p95_yo[wr],
      fillalpha = 0.12, lw = 0, fillcolor = :steelblue, label = "oracle 5–95 %")
plot!(pF, wd, p25_yo[wr]; fillrange = p75_yo[wr],
      fillalpha = 0.28, lw = 0, fillcolor = :steelblue, label = "oracle 25–75 %")

# climatological: 90 % band (light), 50 % band on top (same alpha → appears darker)
plot!(pF, wd, p50_yc[wr]; lw = 2, lc = :darkorange, label = "climatological median")
plot!(pF, wd, p05_yc[wr]; fillrange = p95_yc[wr],
      fillalpha = 0.15, lw = 0, fillcolor = :darkorange, label = "climatological 5–95 %")
plot!(pF, wd, p25_yc[wr]; fillrange = p75_yc[wr],
      fillalpha = 0.35, lw = 0, fillcolor = :darkorange, label = "climatological 25–75 %")

obs_win = [(d, cum_at(days_obs, y_obs, d)) for d in wd if any(abs.(days_obs .- d) .< 4)]
isempty(obs_win) || scatter!(pF, first.(obs_win), last.(obs_win); ms = 5, mc = :black, label = "observed")
savefig(pF, joinpath(figdir, "forecast_$(FORECAST_WKS)wk_day$(fd).png"))
println("Saved: forecast_$(FORECAST_WKS)wk_day$(fd).png")

# =============================================================================
# summary table
# =============================================================================
obs_inc(d1, d2) = cum_at(days_obs, y_obs, d2) - cum_at(days_obs, y_obs, d1)
incr_o = yield_oracle[row(fend), :] .- yield_oracle[row(fd), :]
incr_c = yield_clim[row(fend), :]   .- yield_clim[row(fd), :]
obs_i  = obs_inc(fd, fend)

println("\n=== $(FORECAST_WKS)-week yield forecast  day $fd → $fend ===")
println("                    oracle    climatological")
@printf "  5th  pct :      %5.2f        %5.2f  kg/m²\n" quantile(incr_o,0.05) quantile(incr_c,0.05)
@printf "  median   :      %5.2f        %5.2f  kg/m²\n" quantile(incr_o,0.50) quantile(incr_c,0.50)
@printf "  95th pct :      %5.2f        %5.2f  kg/m²\n" quantile(incr_o,0.95) quantile(incr_c,0.95)
@printf "  observed :      %5.2f\n"                      obs_i
for (label, incr) in [("oracle", incr_o), ("climatological", incr_c)]
    cov = quantile(incr,0.05) <= obs_i <= quantile(incr,0.95) ? "YES ✓" : "NO ✗"
    println("  obs in 90% CI ($label): $cov")
end

println("\nAll figures: ", figdir)
