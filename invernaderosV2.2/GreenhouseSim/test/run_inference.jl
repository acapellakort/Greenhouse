# =============================================================================
# run_inference.jl  --  t-walk MCMC parameter inference on the fast forward map
# =============================================================================
# Run from the GreenhouseSim folder:
#     julia --project=. test/run_inference.jl
#
# Mirrors 2.0's for_inference_V4.jl, but the forward map is the type-stable
# V2.2 one and parameters are selected by name via `update_params` (no globals).
# To infer a different set, just edit `QoI` (any names present in the JSON
# configs, with a prior range added to `QoI_dict`).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Dates, JTwalk, JLD2, Plots

# --- paths -------------------------------------------------------------------
INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))   # -> invernaderos root
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")
obs_path     = joinpath(INV, "invernaderos2.0", "observed_data.csv")

# --- inputs ------------------------------------------------------------------
start_date = DateTime(1998, 7, 11, 0, 0)
end_date   = DateTime(1998, 7, 13, 12, 0)

params   = load_params(climate_json, crop_json)
weather  = load_weather(meteo, start_date, end_date)
controls = load_measured_controls(obs_path)
base     = SimBase(params, weather, controls, datetime2unix(start_date))

days  = [0.0, 1.0]
tspan = (0.0, 86400.0)
teval = 0.0:3600.0:86400.0

# --- observations ------------------------------------------------------------
obs      = DataFrame(CSV.File(obs_path))
tcan_obs = obs.T1
tair_obs = obs.T2
rh_obs   = obs.RH
co2_obs  = obs.C1
error_std = 1.5 .* [0.5, 0.5, 5.0, 0.1]     # Tcan, Tair, RH, CO2

# --- which parameters to infer (edit freely) ---------------------------------
QoI = ["beta2", "gamma3", "gamma4", "nu4"]
QoI_dict = Dict(
    "beta2"  => [0.1,  0.7,   2.0 ],   # [min, true/prior, max]
    "gamma3" => [100.0, 275.0, 500.0],
    "gamma4" => [20.0,  82.0,  200.0],
    "nu4"    => [1e-5,  1e-2,  1e-1 ],
)

chain_length = 250000
burn_in      = 100000

# --- MCMC functions ----------------------------------------------------------
function PriorSupp(x)
    for k in eachindex(x)
        (QoI_dict[QoI[k]][1] < x[k] < QoI_dict[QoI[k]][3]) || return false
    end
    return true
end

function energy(x)
    tcan_p, tair_p, rh_p, co2_p = forward_map(x, QoI, base, days, tspan, teval)
    logL = -sum((tcan_obs .- tcan_p).^2) / error_std[1]^2 -
            sum((tair_obs .- tair_p).^2) / error_std[2]^2 -
            sum((rh_obs   .- rh_p).^2)   / error_std[3]^2 -
            sum((co2_obs  .- co2_p).^2)  / error_std[4]^2
    return -logL
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
savepath = joinpath(postdir, "climate_chain.jld2")
JLD2.@save savepath chain=obj.Output QoI=QoI burn_in=burn_in chain_length=chain_length
println("Climate posterior chain saved to ", savepath)

# --- plots -------------------------------------------------------------------
figdir = joinpath(@__DIR__, "figures")
mkpath(figdir)

# Energy (objective) trace: U = -logL at each accepted state. JTwalk stores it
# in the last column of Output. This is the primary convergence diagnostic --
# after an initial burn-in transient the energy should settle into a stationary
# band with no long-term drift; that plateau is the visual signal of convergence.
energy_trace = obj.Output[:, end]
pE = plot(energy_trace, title = "Energy  U = -logL  (convergence check)",
          xlabel = "iteration", ylabel = "U", label = "")
vline!(pE, [burn_in], lw = 2, lc = :red, ls = :dash, label = "burn-in")
savefig(pE, joinpath(figdir, "energy_trace.png"))

# zoom on the post-burn-in portion (should look like flat noise if converged)
pEz = plot(burn_in:length(energy_trace), energy_trace[burn_in:end],
           title = "Energy (post burn-in)", xlabel = "iteration", ylabel = "U", label = "")
savefig(pEz, joinpath(figdir, "energy_trace_postburnin.png"))

# report basic energy statistics
using Statistics
post = energy_trace[burn_in:end]
println("Energy: min over chain = ", round(minimum(energy_trace), digits = 4))
println("Energy (post burn-in):  mean = ", round(mean(post), digits = 4),
        "  std = ", round(std(post), digits = 4))

# per-parameter chains + posteriors
for (i, name) in enumerate(QoI)
    p1 = plot(obj.Output[:, i], title = "chain: $name", label = "")
    vline!(p1, [burn_in], lw = 2, lc = :red, ls = :dash, label = "burn-in")
    p2 = histogram(obj.Output[burn_in:end, i], title = "posterior: $name", label = "")
    vline!(p2, [QoI_dict[name][2]], lw = 3, lc = :green, label = "prior/true")
    savefig(p1, joinpath(figdir, "$(name)_chain.png"))
    savefig(p2, joinpath(figdir, "$(name)_posterior.png"))
end
println("Done. Figures written to ", figdir)
