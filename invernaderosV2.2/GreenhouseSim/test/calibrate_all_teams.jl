# =============================================================================
# calibrate_all_teams.jl  --  run the growth calibration across all AGC-2018 teams
# =============================================================================
# Run:  julia --project=. test/calibrate_all_teams.jl
#
# Runs the full multi-objective t-walk calibration (weekly yield + leaf number +
# derived LAI) for each grower, at FULL chain length, and produces a cross-team
# comparison of the key posteriors. If the identified parameters (J_max, node_rate,
# veg_sink_max, ...) cluster across independent growers, that is strong evidence the
# model captures real crop physiology rather than fitting one dataset's noise.
#
# Outputs (in test/figures/):
#   teams/<team>_weekly_fit.png, <team>_leaf_fit.png   -- per-team fits
#   teams_params.csv                                   -- mean +/- std per team
#   teams_<param>.png                                  -- cross-team comparison plots

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Dates, JTwalk, JLD2, Plots, Statistics, Printf

# --- config ------------------------------------------------------------------
INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA         = joinpath(INV, "invernaderos2.0", "data")
TEAMS        = ["AiCU", "Croperators", "DeepGreens", "Reference(Growers)", "Sonoma", "iGrow"]
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")
meteo_csv    = joinpath(DATA, "meteo.csv")

chain_length = 200000
burn_in      = 80000

# --- inference specification (same as calibrate_growth.jl) -------------------
QoI_source = ["J_max"]
QoI_growth = ["node_rate", "veg_sink_max", "LAI_max", "set_rate", "Wf_max", "set_start_day"]
QoI_sigma  = ["sigma_y", "sigma_leaf", "sigma_lai"]
QoI        = vcat(QoI_source, QoI_growth, QoI_sigma)
ns, ng     = length(QoI_source), length(QoI_growth)
QoI_dict = Dict(
    "J_max"         => [8e-5,  1.4e-4,  2.8e-4],
    "node_rate"     => [0.02, 0.08, 0.25],
    "veg_sink_max"  => [2.0,  12.0, 40.0],
    "LAI_max"       => [1.5,  2.7,  4.5],
    "set_rate"      => [0.1,  1.0,  5.0],
    "Wf_max"        => [5.0,  16.0, 40.0],
    "set_start_day" => [3.0,  12.0, 30.0],
    "sigma_y"       => [0.02, 0.3,  3.0],
    "sigma_leaf"    => [1.0,  10.0, 40.0],
    "sigma_lai"     => [0.05, 0.5,  2.0],
)

params_base = load_params(climate_json, crop_json)
gp_base     = load_growth_params(growth_json)

# --- helpers -----------------------------------------------------------------
cum_at(days, cum, d) = (d <= days[1] ? float(cum[1]) :
                        d >= days[end] ? float(cum[end]) :
                        let i = searchsortedlast(days, d)
                            days[i] == d ? float(cum[i]) :
                            cum[i] + (d - days[i]) / (days[i+1] - days[i]) * (cum[i+1] - cum[i])
                        end)
weekly_incr(days, cum, edges) = diff([cum_at(days, cum, e) for e in edges])
at_day(v, d) = v[min(d, length(v) - 1) + 1]
safe(name)   = replace(name, r"[^A-Za-z0-9]" => "_")

figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
teamdir = joinpath(figdir, "teams");   mkpath(teamdir)

# --- calibrate one team ------------------------------------------------------
function calibrate_team(TEAM)
    gh_csv   = joinpath(DATA, TEAM, "Greenhouse_climate.csv")
    prod_csv = joinpath(DATA, TEAM, "Production.csv")
    cm_csv   = joinpath(DATA, TEAM, "CropManagement.csv")

    season = load_measured_climate(gh_csv, meteo_csv, params_base)
    days_obs, y_obs = load_production(prod_csv)
    cm_days, cm_Nleaves, cm_LAI = load_cropmanagement(cm_csv)
    edges  = collect(0:7:days_obs[end])
    Δy_obs = weekly_incr(days_obs, y_obs, edges)
    Nwk, Ncm = length(Δy_obs), length(cm_days)

    PriorSupp(z) = all(QoI_dict[QoI[k]][1] < z[k] < QoI_dict[QoI[k]][3] for k in eachindex(z))
    function energy(z)
        p2  = update_params(params_base, QoI_source, z[1:ns])
        gp  = update_params(gp_base,     QoI_growth, z[ns+1:ns+ng])
        σy, σl, σL = z[ns+ng+1], z[ns+ng+2], z[ns+ng+3]
        res = simulate_growth_measured(p2, gp, season; LAI0 = 0.5)
        ry = Δy_obs .- weekly_incr(res.day, res.yield_FW, edges)
        rl = cm_Nleaves .- [at_day(res.node, d) for d in cm_days]
        rL = cm_LAI     .- [at_day(res.LAI,  d) for d in cm_days]
        return Nwk*log(σy) + sum(ry.^2)/(2σy^2) +
               Ncm*log(σl) + sum(rl.^2)/(2σl^2) +
               Ncm*log(σL) + sum(rL.^2)/(2σL^2)
    end

    n   = length(QoI)
    obj = jtwalk(n = n, U = energy, Supp = PriorSupp)
    x0  = [rand()*(QoI_dict[q][3]-QoI_dict[q][1])+QoI_dict[q][1] for q in QoI]
    xp0 = [rand()*(QoI_dict[q][3]-QoI_dict[q][1])+QoI_dict[q][1] for q in QoI]
    Run!(obj, T = chain_length, x0 = x0, xp0 = xp0)

    # save posterior chain --------------------------------------------------
    postdir = joinpath(@__DIR__, "posteriors"); mkpath(postdir)
    JLD2.@save joinpath(postdir, safe(TEAM) * "_growth_chain.jld2") chain=obj.Output QoI=QoI burn_in=burn_in chain_length=chain_length

    pm = [mean(obj.Output[burn_in:end, i]) for i in 1:n]
    ps = [std(obj.Output[burn_in:end, i])  for i in 1:n]

    # per-team fit figures
    p_fit  = update_params(params_base, QoI_source, pm[1:ns])
    gp_fit = update_params(gp_base,     QoI_growth, pm[ns+1:ns+ng])
    res    = simulate_growth_measured(p_fit, gp_fit, season; LAI0 = 0.5)
    s = safe(TEAM)
    pW = plot(edges[2:end], Δy_obs, seriestype=:scatter, ms=4, label="obs",
              xlabel="day", ylabel="kg/m2/wk", title="Weekly yield: $TEAM")
    plot!(pW, edges[2:end], weekly_incr(res.day, res.yield_FW, edges), lw=2, label="model")
    savefig(pW, joinpath(teamdir, "$(s)_weekly_fit.png"))
    pN = plot(cm_days, cm_Nleaves, seriestype=:scatter, ms=4, label="obs",
              xlabel="day", ylabel="leaves/stem", title="Leaf number: $TEAM")
    plot!(pN, res.day, res.node, lw=2, label="model")
    savefig(pN, joinpath(teamdir, "$(s)_leaf_fit.png"))

    yield_final = res.yield_FW[end]; yield_obs = y_obs[end]
    return pm, ps, yield_final, yield_obs
end

# --- run all teams -----------------------------------------------------------
means = Dict{String,Vector{Float64}}()
stds  = Dict{String,Vector{Float64}}()
yfit  = Dict{String,Float64}(); yobs = Dict{String,Float64}()
done  = String[]

for TEAM in TEAMS
    println("\n==================  calibrating ", TEAM, "  ==================")
    try
        @time pm, ps, yf, yo = calibrate_team(TEAM)
        means[TEAM] = pm; stds[TEAM] = ps; yfit[TEAM] = yf; yobs[TEAM] = yo
        push!(done, TEAM)
        @printf("  J_max=%.3e  node_rate=%.4f  veg_sink=%.1f  LAI_max=%.2f  yield %.1f (obs %.1f)\n",
                pm[1], pm[ns+1], pm[ns+2], pm[ns+3], yf, yo)
    catch e
        println("  !! FAILED for ", TEAM, ": ", e)
    end
end

# --- comparison table --------------------------------------------------------
df = DataFrame(team = done)
for (i, q) in enumerate(QoI)
    df[!, q]           = [means[t][i] for t in done]
    df[!, q * "_std"]  = [stds[t][i]  for t in done]
end
df[!, "yield_model"] = [yfit[t] for t in done]
df[!, "yield_obs"]   = [yobs[t] for t in done]
CSV.write(joinpath(figdir, "teams_params.csv"), df)

# --- cross-team comparison plots ---------------------------------------------
compare = ["J_max", "node_rate", "veg_sink_max", "LAI_max", "set_start_day"]
idx(q)  = findfirst(==(q), QoI)
for q in compare
    i = idx(q)
    m = [means[t][i] for t in done]
    s = [stds[t][i]  for t in done]
    p = scatter(1:length(done), m, yerror = s, ms = 6, legend = false,
                xticks = (1:length(done), done), xrotation = 30,
                ylabel = q, title = "$q across teams (mean +/- std)")
    savefig(p, joinpath(figdir, "teams_$(q).png"))
end

println("\nAll done. Calibrated: ", join(done, ", "))
println("Table: ", joinpath(figdir, "teams_params.csv"))
println("Per-param comparison plots: teams_<param>.png ;  per-team fits in teams/")
