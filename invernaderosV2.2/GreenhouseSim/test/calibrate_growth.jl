# =============================================================================
# calibrate_growth.jl  --  MULTI-OBJECTIVE cucumber growth calibration
# =============================================================================
# Run:  julia --project=. test/calibrate_growth.jl
#
# Photosynthesis is driven by MEASURED greenhouse climate (fast, no ODE). We fit
# THREE observables jointly, each with its own inferred sigma:
#   1. weekly incremental yield        (Production.csv, Total_Prod_cum)
#   2. cumulative leaf/node number     (CropManagement.csv, N_leaves)  -> development
#   3. derived standing LAI            (from N_leaves)                  -> vegetative allocation
#
# FIXED PARAMETERS (not inferred):
#   J_max         = 2.0e-4 mol e⁻ m⁻² s⁻¹   AiCU is light-limited throughout;
#                    J_max is unidentifiable from this dataset (flat posterior at 200k steps).
#   Wf_max        = 16.0 g DW / fruit  400 g FW (commercial high-wire standard)
#                    × DMC_fruit 0.04 = 16 g DW.  Encodes the HARVEST DECISION;
#                    team-specific: Wf_max = FW_harvest × DMC_fruit.
#   set_start_day = 15 days             grower observation (first fruit seen Aug 29 = day 15).
#                    Sweep over {12..16} confirms day 15 minimises full-season yield RSS
#                    and matches the posterior mean (15.05 d) from the 8-param run.
#
# INFERRED (7 parameters):
#   node_rate, veg_sink_max, LAI_max, set_rate  (growth)
#   sigma_y, sigma_leaf, sigma_lai              (noise)
#
# ASSUMPTIONS in the LAI derivation:
#   leaf_area = 0.05 m^2/leaf, stem_density = 2.5 stems/m^2, N_window = 20 leaves.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Dates, JTwalk, JLD2, Plots, Statistics

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

# --- load model + data -------------------------------------------------------
params  = load_params(climate_json, crop_json)
gp_base = load_growth_params(growth_json)   # Wf_max = 16.0 fixed in JSON
gp_base = update_params(gp_base, ["set_start_day"], [15.0])  # fixed from grower observation + RSS sweep (diagnose_set_start_day.jl)

println("Loading measured climate (", TEAM, ") ...")
season = load_measured_climate(gh_csv, meteo_csv, params)
println("  season length = ", season.ndays, " days")

# observable 1: weekly incremental yield
days_obs, y_obs = load_production(prod_csv)
cum_at(days, cum, d) = (d <= days[1] ? float(cum[1]) :
                        d >= days[end] ? float(cum[end]) :
                        let i = searchsortedlast(days, d)
                            days[i] == d ? float(cum[i]) :
                            cum[i] + (d - days[i]) / (days[i+1] - days[i]) * (cum[i+1] - cum[i])
                        end)
weekly_incr(days, cum, edges) = diff([cum_at(days, cum, e) for e in edges])
edges  = collect(0:7:days_obs[end])
Δy_obs = weekly_incr(days_obs, y_obs, edges)
Nwk    = length(Δy_obs)

# observables 2 & 3: leaf number and derived LAI
cm_days, cm_Nleaves, cm_LAI = load_cropmanagement(cm_csv)
Ncm = length(cm_days)
at_day(v, d) = v[min(d, length(v) - 1) + 1]
println("  ", Nwk, " weekly yield pts, ", Ncm, " crop-management pts")
println("  observed final yield = ", round(y_obs[end], digits = 1), " kg/m2; ",
        "final leaves ~ ", round(cm_Nleaves[end]), "; LAI plateau ~ ",
        round(cm_LAI[end], digits = 2))


# --- parameters to infer -----------------------------------------------------
# Wf_max FIXED: 400 g FW × DMC 0.04 = 16 g DW (already set in cucumber_growth.json)
# J_max  FIXED: 2.0e-4 mol e⁻ m⁻² s⁻¹ (unidentifiable under AiCU light-limited conditions)
QoI_growth = ["node_rate", "veg_sink_max", "LAI_max", "set_rate"]
QoI_sigma  = ["sigma_y", "sigma_leaf", "sigma_lai"]
QoI        = vcat(QoI_growth, QoI_sigma)
ng         = length(QoI_growth)   # 4
n          = length(QoI)          # 7

QoI_dict = Dict(
    "node_rate"     => [0.02, 0.08,  0.25],   # nodes per degree-day
    "veg_sink_max"  => [2.0,  12.0,  40.0],   # g DM m-2 d-1
    "LAI_max"       => [1.5,  2.7,   4.5],    # canopy LAI cap
    "set_rate"      => [0.1,  1.0,   5.0],    # fruits m-2 d-1
    "set_start_day" => [3.0,  12.0,  30.0],   # fruit-set lag (days)
    "sigma_y"       => [0.02, 0.3,   3.0],    # kg m-2 wk-1
    "sigma_leaf"    => [1.0,  10.0,  40.0],   # leaves
    "sigma_lai"     => [0.05, 0.5,   2.0],    # LAI units
)
chain_length = 200000
burn_in      = 80000

# --- MCMC functions ----------------------------------------------------------
function PriorSupp(z)
    for k in eachindex(z)
        (QoI_dict[QoI[k]][1] < z[k] < QoI_dict[QoI[k]][3]) || return false
    end
    return true
end

function energy(z)
    gp  = update_params(gp_base, QoI_growth, z[1:ng])
    σy, σl, σL = z[ng+1], z[ng+2], z[ng+3]
    res = simulate_growth_measured(params, gp, season; LAI0 = 0.5)

    ry = Δy_obs .- weekly_incr(res.day, res.yield_FW, edges)
    rl = cm_Nleaves .- [at_day(res.node, d) for d in cm_days]
    rL = cm_LAI     .- [at_day(res.LAI,  d) for d in cm_days]

    return Nwk * log(σy) + sum(ry .^ 2) / (2σy^2) +
           Ncm * log(σl) + sum(rl .^ 2) / (2σl^2) +
           Ncm * log(σL) + sum(rL .^ 2) / (2σL^2)
end

# --- run ---------------------------------------------------------------------
obj = jtwalk(n = n, U = energy, Supp = PriorSupp)
x0  = [rand() * (QoI_dict[q][3] - QoI_dict[q][1]) + QoI_dict[q][1] for q in QoI]
xp0 = [rand() * (QoI_dict[q][3] - QoI_dict[q][1]) + QoI_dict[q][1] for q in QoI]

println("Inferring: ", QoI)
println("  Fixed: Wf_max = $(gp_base.Wf_max) g DW/fruit  (400 g FW × DMC 0.04)")
println("  Fixed: set_start_day = 15 d  (grower obs + RSS sweep)")
@time Run!(obj, T = chain_length, x0 = x0, xp0 = xp0)

# --- save chain --------------------------------------------------------------
postdir  = joinpath(@__DIR__, "posteriors"); mkpath(postdir)
savepath = joinpath(postdir, TEAM * "_growth_chain.jld2")
JLD2.@save savepath chain=obj.Output QoI=QoI burn_in=burn_in chain_length=chain_length
println("Posterior chain saved to ", savepath)

# --- posterior summary -------------------------------------------------------
figdir    = joinpath(@__DIR__, "figures"); mkpath(figdir)
post_samp = obj.Output[burn_in:end, :]
post_mean = vec(mean(post_samp, dims = 1))
post_std  = vec(std(post_samp,  dims = 1))
println("\n--- posterior (mean ± std) ---")
for (i, name) in enumerate(QoI)
    println("  ", rpad(name, 14), " = ", round(post_mean[i], digits = 4),
            "  ±  ", round(post_std[i], digits = 4))
end
println("  (fixed) Wf_max       = ", gp_base.Wf_max, " g DW/fruit")
println("  (fixed) set_start_day = 15 days  (grower observation + RSS sweep)")


pE = plot(obj.Output[:, end], title = "Energy U = -logL",
          xlabel = "iteration", ylabel = "U", label = "")
vline!(pE, [burn_in], lw = 2, lc = :red, ls = :dash, label = "burn-in")
savefig(pE, joinpath(figdir, "growthcal_energy.png"))

for (i, name) in enumerate(QoI)
    p1 = plot(obj.Output[:, i], title = "chain: $name", label = "")
    vline!(p1, [burn_in], lw = 2, lc = :red, ls = :dash, label = "burn-in")
    p2 = histogram(post_samp[:, i], title = "posterior: $name", label = "")
    vline!(p2, [post_mean[i]], lw = 3, lc = :green, label = "mean")
    savefig(p1, joinpath(figdir, "growthcal_$(name)_chain.png"))
    savefig(p2, joinpath(figdir, "growthcal_$(name)_posterior.png"))
end

# --- fit at posterior mean ---------------------------------------------------
gp_fit  = update_params(gp_base, QoI_growth, post_mean[1:ng])
res_fit = simulate_growth_measured(params, gp_fit, season; LAI0 = 0.5)

pW = plot(edges[2:end], Δy_obs, seriestype = :scatter, ms = 4, label = "observed",
          xlabel = "day", ylabel = "kg FW/m2/wk", title = "Weekly yield ($TEAM)")
plot!(pW, edges[2:end], weekly_incr(res_fit.day, res_fit.yield_FW, edges), lw = 2, label = "model")
savefig(pW, joinpath(figdir, "growthcal_weekly_fit.png"))

pC = plot(days_obs, y_obs, lw = 2, label = "observed",
          xlabel = "day", ylabel = "kg FW/m2", title = "Cumulative yield ($TEAM)")
plot!(pC, res_fit.day, res_fit.yield_FW, lw = 2, ls = :dash, label = "model")
savefig(pC, joinpath(figdir, "growthcal_cumulative_fit.png"))

pN = plot(cm_days, cm_Nleaves, seriestype = :scatter, ms = 4, label = "observed",
          xlabel = "day", ylabel = "leaves / stem", title = "Leaf number ($TEAM)")
plot!(pN, res_fit.day, res_fit.node, lw = 2, label = "model")
savefig(pN, joinpath(figdir, "growthcal_leaf_fit.png"))

pL = plot(cm_days, cm_LAI, seriestype = :scatter, ms = 4, label = "observed (derived)",
          xlabel = "day", ylabel = "LAI", title = "LAI ($TEAM)")
plot!(pL, res_fit.day, res_fit.LAI, lw = 2, label = "model")
savefig(pL, joinpath(figdir, "growthcal_lai_fit.png"))

println("\nDone. Figures in ", figdir)
println("Model final yield = ", round(res_fit.yield_FW[end], digits = 2),
        " kg/m2  (obs ", round(y_obs[end], digits = 2), ")")
println("set_rate posterior: ", round(post_mean[findfirst(==("set_rate"), QoI)], digits = 3),
        " ± ", round(post_std[findfirst(==("set_rate"), QoI)], digits = 3), " fruits m⁻² d⁻¹")
