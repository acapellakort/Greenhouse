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
# Rationale: yield ALONE only constrains the net source, leaving the partitioning
# parameters (veg_sink_max, set_rate, Wf_max) unidentified and forcing the light-
# extinction k_ext to an unphysical value to absorb a source bias. Adding the leaf
# / LAI observables constrains development and vegetative growth independently, so
# k_ext can be FIXED at a physical 0.75 and the partitioning becomes identifiable.
#
# ASSUMPTIONS in the LAI derivation (no measured LAI / plant density in the data):
#   leaf_area = 0.05 m^2/leaf, stem_density = 2.5 stems/m^2, N_window = 20 leaves.
#   Change these at the load_cropmanagement call if the real values are known.

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
gp_base = load_growth_params(growth_json)          # k_ext=0.75, SLA=0.03 fixed here

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
cm_days, cm_Nleaves, cm_LAI = load_cropmanagement(cm_csv)   # literature defaults; see header
Ncm = length(cm_days)
at_day(v, d) = v[min(d, length(v) - 1) + 1]
println("  ", Nwk, " weekly yield pts, ", Ncm, " crop-management pts")
println("  observed final yield = ", round(y_obs[end], digits = 1), " kg/m2; ",
        "final leaves ~ ", round(cm_Nleaves[end]), "; LAI plateau ~ ", round(cm_LAI[end], digits = 2))

# --- parameters to infer (k_ext & SLA are FIXED in the JSON) ------------------
# k_ext is FIXED at the physical 0.75 in the JSON. We infer only J_max: with AiCU's
# heavy CO2 dosing (633-2000 ppm) the crop is LIGHT-limited throughout, so the FvCB
# min(Ac, Aj, Ap) is set by the electron-transport limit Aj -> only J_max binds.
# V_cmax25 is NOT identifiable from this dataset (flat posterior when inferred), so
# it is fixed at its literature value in the crop JSON (does not affect the fit
# while it is non-binding). It would only become identifiable under Rubisco-limiting
# data (low CO2, high light).
QoI_source = ["J_max"]                                                # in `params`
QoI_growth = ["node_rate", "veg_sink_max", "LAI_max", "set_rate", "Wf_max", "set_start_day"]  # in `gp`
QoI_sigma  = ["sigma_y", "sigma_leaf", "sigma_lai"]
QoI        = vcat(QoI_source, QoI_growth, QoI_sigma)
ns         = length(QoI_source)
ng         = length(QoI_growth)
QoI_dict = Dict(
    "J_max"         => [8e-5,  1.4e-4,  2.8e-4],  # mol e- m-2 s-1  (was 2.0e-4, generic)
    "node_rate"     => [0.02, 0.08, 0.25],   # nodes per degree-day
    "veg_sink_max"  => [2.0,  12.0, 40.0],   # g DM/m2/day
    "LAI_max"       => [1.5,  2.7,  4.5],     # canopy LAI cap
    "set_rate"      => [0.1,  1.0,  5.0],     # fruits/m2/day
    "Wf_max"        => [5.0,  16.0, 40.0],    # g DW/fruit
    "set_start_day" => [3.0,  12.0, 30.0],    # fruit-set lag (days)
    "sigma_y"       => [0.02, 0.3,  3.0],     # kg/m2/week
    "sigma_leaf"    => [1.0,  10.0, 40.0],    # leaves
    "sigma_lai"     => [0.05, 0.5,  2.0],     # LAI units
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
    p2  = update_params(params,  QoI_source, z[1:ns])              # V_cmax25, J_max
    gp  = update_params(gp_base, QoI_growth, z[ns+1:ns+ng])
    σy, σl, σL = z[ns+ng+1], z[ns+ng+2], z[ns+ng+3]
    res = simulate_growth_measured(p2, gp, season; LAI0 = 0.5)

    ry = Δy_obs .- weekly_incr(res.day, res.yield_FW, edges)
    rl = cm_Nleaves .- [at_day(res.node, d) for d in cm_days]
    rL = cm_LAI     .- [at_day(res.LAI,  d) for d in cm_days]

    return Nwk  * log(σy) + sum(ry .^ 2) / (2σy^2) +
           Ncm  * log(σl) + sum(rl .^ 2) / (2σl^2) +
           Ncm  * log(σL) + sum(rL .^ 2) / (2σL^2)
end

# --- run ---------------------------------------------------------------------
n   = length(QoI)
obj = jtwalk(n = n, U = energy, Supp = PriorSupp)
x0  = [rand() * (QoI_dict[q][3] - QoI_dict[q][1]) + QoI_dict[q][1] for q in QoI]
xp0 = [rand() * (QoI_dict[q][3] - QoI_dict[q][1]) + QoI_dict[q][1] for q in QoI]

println("Inferring: ", QoI)
@time Run!(obj, T = chain_length, x0 = x0, xp0 = xp0)

# --- save posterior chain (JLD2) ---------------------------------------------
postdir = joinpath(@__DIR__, "posteriors"); mkpath(postdir)
savepath = joinpath(postdir, TEAM * "_growth_chain.jld2")
JLD2.@save savepath chain=obj.Output QoI=QoI burn_in=burn_in chain_length=chain_length
println("Posterior chain saved to ", savepath)

# --- posterior summary -------------------------------------------------------
figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
post_mean = [mean(obj.Output[burn_in:end, i]) for i in 1:n]
post_std  = [std(obj.Output[burn_in:end, i])  for i in 1:n]
println("\n--- posterior (mean +/- std) ---")
for (i, name) in enumerate(QoI)
    println("  ", rpad(name, 14), " = ", round(post_mean[i], digits = 4),
            "  +/- ", round(post_std[i], digits = 4))
end

energy_trace = obj.Output[:, end]
pE = plot(energy_trace, title = "Energy U = -logL (multi-objective)",
          xlabel = "iteration", ylabel = "U", label = "")
vline!(pE, [burn_in], lw = 2, lc = :red, ls = :dash, label = "burn-in")
savefig(pE, joinpath(figdir, "growthcal_energy.png"))

for (i, name) in enumerate(QoI)
    p1 = plot(obj.Output[:, i], title = "chain: $name", label = "")
    vline!(p1, [burn_in], lw = 2, lc = :red, ls = :dash, label = "burn-in")
    p2 = histogram(obj.Output[burn_in:end, i], title = "posterior: $name", label = "")
    vline!(p2, [post_mean[i]], lw = 3, lc = :green, label = "mean")
    savefig(p1, joinpath(figdir, "growthcal_$(name)_chain.png"))
    savefig(p2, joinpath(figdir, "growthcal_$(name)_posterior.png"))
end

# --- fits at posterior mean --------------------------------------------------
p_fit   = update_params(params,  QoI_source, post_mean[1:ns])
gp_fit  = update_params(gp_base, QoI_growth, post_mean[ns+1:ns+ng])
res_fit = simulate_growth_measured(p_fit, gp_fit, season; LAI0 = 0.5)

pW = plot(edges[2:end], Δy_obs, seriestype = :scatter, ms = 4, label = "observed",
          xlabel = "day", ylabel = "kg FW/m2/wk", title = "Weekly yield ($TEAM)")
plot!(pW, edges[2:end], weekly_incr(res_fit.day, res_fit.yield_FW, edges), lw = 2, label = "model")
savefig(pW, joinpath(figdir, "growthcal_weekly_fit.png"))

pC = plot(days_obs, y_obs, lw = 2, label = "observed", xlabel = "day",
          ylabel = "kg FW/m2", title = "Cumulative yield ($TEAM)")
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
        " kg/m2 (obs ", round(y_obs[end], digits = 2), ");  final leaves = ",
        round(res_fit.node[end]), " (obs ", round(cm_Nleaves[end]), ")")
