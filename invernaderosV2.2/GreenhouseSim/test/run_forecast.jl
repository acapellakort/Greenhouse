# =============================================================================
# run_forecast.jl  --  Probabilistic 3–5 week production forecast
# =============================================================================
# Run from the GreenhouseSim folder:
#     julia --project=. test/run_forecast.jl
#
# Uses the t-walk posterior saved by calibrate_growth.jl to propagate parameter
# uncertainty forward from a given crop day (FORECAST_DAY) for HORIZON days.
# Weather uncertainty is captured by replaying the forecast window under each of
# WEATHER_YEARS historical weather series.  Outputs fan charts + marginal
# distributions saved to test/figures/.
#
# Prerequisites:
#   julia --project=. test/calibrate_growth.jl
#   (creates test/posteriors/<TEAM>_growth_chain.jld2)
#
# Ground rules (project notes):
#   - V2.2 engine frozen; physics-flag overrides stay default-off.
#   - Forecast weather drawn from historical years ≠ calibration year (2018).
#   - Caveats: the forecast uses big-leaf photosynthesis (not cohort); cohort
#     twin should be re-calibrated separately before using cohort forecast.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Dates, JLD2, Plots, Statistics
import Base.deepcopy

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
TEAM          = "AiCU"
FORECAST_DAY  = 45          # crop day from which to forecast (0-indexed)
HORIZON       = 35          # forecast horizon in days (5 weeks)
N_SAMPLES     = 200         # posterior samples to use (thinned)
# Weather years for forecast: 1989–2017 (excludes calibration year 2018 and
# corrupt 2023; 1988 is partial starting July so skip it).
WEATHER_YEARS = collect(1989:2017)   # 29 historical scenarios
# AGC season start month/day (Aug 14 in 2018):
SEASON_MONTH  = 8
SEASON_DAY    = 14

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA         = joinpath(INV, "invernaderos2.0", "data")
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")
meteo_csv    = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")
gh_csv       = joinpath(DATA, TEAM, "Greenhouse_climate.csv")
obs_path     = normpath(joinpath(INV, "invernaderos2.0", "observed_data.csv"))
prod_csv     = joinpath(DATA, TEAM, "Production.csv")
cm_csv       = joinpath(DATA, TEAM, "CropManagement.csv")
chain_file   = joinpath(@__DIR__, "posteriors", TEAM * "_growth_chain.jld2")

# ---------------------------------------------------------------------------
# Load posterior chain
# ---------------------------------------------------------------------------
isfile(chain_file) || error("""
  Posterior chain not found at:
    $chain_file
  Please run first:
    julia --project=. test/calibrate_growth.jl
""")

JLD2.@load chain_file chain QoI burn_in chain_length
post = chain[burn_in+1:end, :]         # rows after burn-in
thin = max(1, size(post, 1) ÷ N_SAMPLES)
post = post[1:thin:end, :]             # thinned
N    = size(post, 1)
println("Posterior samples: ", N, "  (thinned from ", size(chain, 1) - burn_in, " post-burn-in)")

ns = 1          # J_max is the only source-side QoI
ng = 6          # growth QoIs (see calibrate_growth.jl)
# QoI = ["J_max", "node_rate", "veg_sink_max", "LAI_max", "set_rate", "Wf_max",
#         "set_start_day", "sigma_y", "sigma_leaf", "sigma_lai"]

# ---------------------------------------------------------------------------
# Base params + measured season (for phase-1: 0..FORECAST_DAY)
# ---------------------------------------------------------------------------
params_base = load_params(climate_json, crop_json)
gp_base     = load_growth_params(growth_json)

println("Loading measured calibration climate (", TEAM, ") ...")
season = load_measured_climate(gh_csv, joinpath(DATA, "meteo.csv"), params_base)
println("  season length = ", season.ndays, " days")
FORECAST_DAY <= season.ndays || error("FORECAST_DAY=$FORECAST_DAY exceeds season length $(season.ndays)")

days_obs, y_obs = load_production(prod_csv)
cm_days, cm_Nleaves, cm_LAI = load_cropmanagement(cm_csv)

# ---------------------------------------------------------------------------
# Phase 1: run 0..FORECAST_DAY under measured climate for each posterior sample
# ---------------------------------------------------------------------------
println("Running phase-1 (measured climate, 0..", FORECAST_DAY, " days) for ", N, " samples ...")

# helper: run phase-1 under measured climate, return state + history
function phase1_run(p, gp, season, D)
    s = init_growth_state(gp; LAI0 = 0.5)
    yield_hist = Float64[]
    LAI_hist   = Float64[]
    node_hist  = Float64[]
    for d in 0:(D - 1)
        M   = season.day_samples[d + 1]
        Pgd = daily_assimilation_measured(M, s.LAI, p, season.dt, gp.k_ext)
        grow!(s, Pgd, GreenhouseSim.daily_mean_T(M), gp, d)
        push!(yield_hist, s.yield_FW)
        push!(LAI_hist,   s.LAI)
        push!(node_hist,  s.node)
    end
    return s, yield_hist, LAI_hist, node_hist
end

states_at_D  = Vector{GrowthState}(undef, N)
yield_phase1 = Matrix{Float64}(undef, N, FORECAST_DAY)
LAI_phase1   = Matrix{Float64}(undef, N, FORECAST_DAY)
node_phase1  = Matrix{Float64}(undef, N, FORECAST_DAY)

for i in 1:N
    p  = update_params(params_base,  ["J_max"], post[i, 1:ns])
    gp = update_params(gp_base, ["node_rate","veg_sink_max","LAI_max","set_rate","Wf_max","set_start_day"],
                       post[i, ns+1:ns+ng])
    s, yh, Lh, nh = phase1_run(p, gp, season, FORECAST_DAY)
    states_at_D[i] = s
    yield_phase1[i, :] = yh
    LAI_phase1[i, :] = Lh
    node_phase1[i, :] = nh
    i % 50 == 0 && println("  phase-1 sample ", i, "/", N)
end
println("Phase-1 done.")

# ---------------------------------------------------------------------------
# Phase 2: forecast forward HORIZON days under each historical weather year
# ---------------------------------------------------------------------------
# AGC season starts Aug-14 each year; we take the same calendar window.
println("Loading meteo (", length(WEATHER_YEARS), " weather scenarios) ...")
meteo_df = GreenhouseSim.prepare_meteo(meteo_csv)

# Measured controls (AiCU); extrapolate beyond the calibration period
controls = load_measured_controls(obs_path)

# Forward-integrate one trajectory: state s0, params (p,gp), weather, starting
# at day D_start of the season (for `grow!` day index continuity), for H days.
function phase2_run(s0::GrowthState, p, gp, weather, controls, season_start_unix, D_start, H)
    s = deepcopy(s0)
    yields = Float64[]
    LAIs   = Float64[]
    nodes  = Float64[]
    u = [18 + 273.15, 23 + 273.15, 575.0, 1200.0, 0.0, 0.0, 18 + 273.15]
    tspan = (0.0, 86400.0)
    teval = 0.0:600.0:86400.0
    for h in 0:(H - 1)
        d        = D_start + h
        wt0      = season_start_unix + 86400.0 * d
        ct0      = 86400.0 * d
        u[6]     = 0.0
        Pg, Tmean, u = day_assimilation(p, weather, controls, s.LAI, wt0, ct0, u, tspan, teval)
        grow!(s, Pg, Tmean, gp, d)
        push!(yields, s.yield_FW)
        push!(LAIs,   s.LAI)
        push!(nodes,  s.node)
    end
    return yields, LAIs, nodes
end

# Collect trajectories: [sample, weather_year, horizon_day]
nW = length(WEATHER_YEARS)
forecast_yield = Array{Float64}(undef, N, nW, HORIZON)
forecast_LAI   = Array{Float64}(undef, N, nW, HORIZON)

println("Running phase-2 forecast (", N, " samples × ", nW, " weather years × ", HORIZON, " days) ...")
for (wi, yr) in enumerate(WEATHER_YEARS)
    # Construct weather for this year's Aug-14 .. Aug-14 + FORECAST_DAY + HORIZON + 5
    t0 = DateTime(yr, SEASON_MONTH, SEASON_DAY, 0, 0)
    t1 = t0 + Day(FORECAST_DAY + HORIZON + 5)
    local weather
    try
        weather = weather_from_df(meteo_df, t0, t1)
    catch e
        println("  Skipping year $yr (window missing): $e")
        forecast_yield[:, wi, :] .= NaN
        forecast_LAI[:, wi, :]   .= NaN
        continue
    end
    season_start_unix = datetime2unix(t0)

    for i in 1:N
        p  = update_params(params_base, ["J_max"], post[i, 1:ns])
        gp = update_params(gp_base, ["node_rate","veg_sink_max","LAI_max","set_rate","Wf_max","set_start_day"],
                           post[i, ns+1:ns+ng])
        yf, Lf, _ = phase2_run(states_at_D[i], p, gp, weather, controls,
                                season_start_unix, FORECAST_DAY, HORIZON)
        # yields here are cumulative from season start including phase-1 yield
        forecast_yield[i, wi, :] = yf  # cumulative from day-0 (deepcopy carries phase-1 yield)
        forecast_LAI[i, wi, :]   = Lf
    end
    wi % 5 == 0 && println("  weather year $yr ($wi/$nW) done")
end
println("Phase-2 done.")

# ---------------------------------------------------------------------------
# Summarise: quantiles across (sample × weather_year) ensemble
# ---------------------------------------------------------------------------
# Flatten sample and weather dimensions
flat_yield = reshape(forecast_yield, N * nW, HORIZON)   # rows: trajectories
flat_LAI   = reshape(forecast_LAI,   N * nW, HORIZON)

# Filter out any NaN trajectories
valid = .!any(isnan, flat_yield, dims=2)[:, 1]
flat_yield = flat_yield[valid, :]
flat_LAI   = flat_LAI[valid, :]
println("Valid trajectories: ", sum(valid), " / ", N * nW)

qs = [0.05, 0.25, 0.50, 0.75, 0.95]
quant_yield = mapslices(v -> quantile(v, qs), flat_yield; dims=1)  # (5, HORIZON)
quant_LAI   = mapslices(v -> quantile(v, qs), flat_LAI;   dims=1)

forecast_days = FORECAST_DAY .+ (1:HORIZON)

# Phase-1 quantiles (across posterior samples only, no weather uncertainty)
p1_yield_q = mapslices(v -> quantile(v, qs), yield_phase1; dims=1)
p1_LAI_q   = mapslices(v -> quantile(v, qs), LAI_phase1;   dims=1)
phase1_days = 1:FORECAST_DAY

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
clrs = [:blue, :orange]

function fan_plot!(plt, days, Q; alpha_lo=0.15, alpha_hi=0.30, color=:steelblue, label_median="median")
    plot!(plt, days, Q[3, :], lw=2, lc=color, label=label_median)
    plot!(plt, days, Q[2, :], fillrange=Q[4, :], fillalpha=alpha_hi, lc=:transparent,
          fillcolor=color, label="50% CI")
    plot!(plt, days, Q[1, :], fillrange=Q[5, :], fillalpha=alpha_lo, lc=:transparent,
          fillcolor=color, label="90% CI")
end

# --- Yield fan chart ---
pY = plot(title = "Cumulative yield forecast — $TEAM (start day $FORECAST_DAY)",
          xlabel = "Crop day", ylabel = "Yield [kg FW / m²]", legend = :topleft)
# Phase-1 (measured climate, parameter uncertainty only)
fan_plot!(pY, phase1_days, p1_yield_q; color=:steelblue, label_median="hist. median")
# Phase-2 (forecast, parameter + weather uncertainty)
fan_plot!(pY, forecast_days, quant_yield; color=:firebrick, label_median="forecast median")
# Forecast start marker
vline!(pY, [FORECAST_DAY], lw=1.5, lc=:black, ls=:dash, label="forecast start")
# Observations
obs_in  = days_obs[days_obs .<= FORECAST_DAY + HORIZON .+ 3]
y_in    = y_obs[days_obs .<= FORECAST_DAY + HORIZON .+ 3]
scatter!(pY, obs_in, y_in, ms=5, mc=:black, label="observed")
savefig(pY, joinpath(figdir, "forecast_$(TEAM)_yield.png"))

# --- LAI fan chart ---
pL = plot(title = "LAI forecast — $TEAM (start day $FORECAST_DAY)",
          xlabel = "Crop day", ylabel = "LAI", legend = :topleft)
fan_plot!(pL, phase1_days, p1_LAI_q;  color=:steelblue, label_median="hist. median")
fan_plot!(pL, forecast_days, quant_LAI; color=:firebrick, label_median="forecast median")
vline!(pL, [FORECAST_DAY], lw=1.5, lc=:black, ls=:dash, label="forecast start")
obs_LAI_in = cm_LAI[cm_days .<= FORECAST_DAY + HORIZON .+ 3]
days_LAI_in = cm_days[cm_days .<= FORECAST_DAY + HORIZON .+ 3]
scatter!(pL, days_LAI_in, obs_LAI_in, ms=5, mc=:black, label="observed (derived)")
savefig(pL, joinpath(figdir, "forecast_$(TEAM)_LAI.png"))

# --- Marginal yield distribution at horizon ---
final_yields = flat_yield[:, end]
pD = histogram(final_yields, nbins=30, normalize=:pdf, xlabel="Yield [kg FW/m²]",
               ylabel="Density", title="Forecast yield distribution at day $(FORECAST_DAY + HORIZON) ($TEAM)",
               label="")
vline!(pD, [quantile(final_yields, 0.05)], lw=2, lc=:red, label="5th pct")
vline!(pD, [quantile(final_yields, 0.50)], lw=2, lc=:green, label="median")
vline!(pD, [quantile(final_yields, 0.95)], lw=2, lc=:red, ls=:dash, label="95th pct")
if !isempty(y_obs)
    vline!(pD, [y_obs[end]], lw=2, lc=:black, label="final obs")
end
savefig(pD, joinpath(figdir, "forecast_$(TEAM)_marginal.png"))

# Combined figure
savefig(plot(pY, pL, pD, layout=(1, 3), size=(1400, 450)),
        joinpath(figdir, "forecast_$(TEAM)_combined.png"))

# ---------------------------------------------------------------------------
# Numerical summary
# ---------------------------------------------------------------------------
println("\n--- Forecast summary ($TEAM, day $FORECAST_DAY → $(FORECAST_DAY + HORIZON)) ---")
println("Yield at day $(FORECAST_DAY + HORIZON):")
for (q, v) in zip(["5th", "25th", "50th", "75th", "95th"], quant_yield[:, end])
    println("  ", rpad(q, 5), " pct : ", round(v, digits=2), " kg/m²")
end
println("LAI at day $(FORECAST_DAY + HORIZON):")
for (q, v) in zip(["5th", "25th", "50th", "75th", "95th"], quant_LAI[:, end])
    println("  ", rpad(q, 5), " pct : ", round(v, digits=3))
end
println("\nFigures written to ", figdir)
println("\nCaveats:")
println("  - Big-leaf photosynthesis; cohort twin gives +8% yield and should be")
println("    re-calibrated separately (see notes/14 and project F2).")
println("  - Phase-1 state uncertainty propagated; phase-2 weather drawn from")
println("    historical years 1989-2017 (not operational NWP).")
println("  - Controls held at AiCU measured values throughout; forecast is for")
println("    AiCU management strategy only.")
