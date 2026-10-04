using GreenhouseSim
using CSV
using DataFrames
using Dates
using DifferentialEquations
using JTwalk
using Plots

# 1. --- Basic Configuration ---
start_date = DateTime(1998, 7, 11, 0, 0)
No_days    = 1  # Number of days to simulate (2-1 in your original)
days       = 0:No_days
t_span     = (0.0, 24 * 3600.0)
t_eval     = (0.0:3600.0:24*3600.0) # Reporting every hour

# 2. --- File Paths ---
base_dir      = pwd()
data_path     = joinpath(base_dir, "data", "observed_data.csv")
climate_json  = joinpath(base_dir, "configfiles", "constants_climate.json")
crop_json     = joinpath(base_dir, "configfiles", "constants_cropphoto.json")
controls_json = joinpath(base_dir, "configfiles", "control_instructions.json")

# 3. --- Initialize Data Containers ---
# Load Parameters (Uses the unit_mapping defined in your package)
# Note: We merge climate and crop parameters into one ModelParams object
p_climate = load_model_params(climate_json, GreenhouseSim.unit_mapping)
p_crop    = load_model_params(crop_json, GreenhouseSim.unit_mapping)
params    = GreenhouseSim.merge_params(p_climate, p_crop)

# Load Controls
cs_json = load_controls(controls_json)
cs_csv  = load_control_interpolations(data_path)

# Load Weather (Assuming you have a loader in your package)
# If not, you can keep your original weather interpolation logic here
weather_dict = GreenhouseSim.load_weather_data(joinpath(base_dir, "utils", "weather_data.csv"))

# 4. --- Create Simulation Context ---
# This object replaces all your previous 'global' variables
base_ctx = SimulationContext(
    params = params,
    controls_json = cs_json,
    controls_csv = cs_csv,
    weather = weather_dict,
    start_time_unix = datetime2unix(start_date),
    start_date = start_date
)

# 5. --- Observed Data Setup ---
observed_data = DataFrame(CSV.File(data_path))
# Standard deviations for Likelihood
error_std = [0.75, 0.75, 7.5, 0.15] # Tcan, Tair, RH, CO2

# Extract target vectors
tcan_obs = observed_data.T1 
tair_obs = observed_data.T2 
rh_obs   = observed_data.RH  
co2_obs  = observed_data.C1 

# 6. --- Inference Setup ---
QoI = ["beta2", "gamma3", "gamma4", "nu4"]
QoI_dict = Dict(
    "beta2"     => [0.1, 0.7, 2.0],
    "gamma3"    => [100.0, 275.0, 500.0],
    "gamma4"    => [20.0, 82.0, 200.0],
    "nu4"       => [1e-5, 1e-2, 1e-1]
)

chain_length = 2500
burn_in      = 2000

# 7. --- MCMC Functions ---

function PriorSupp(x)
    for i in 1::Integer:length(x)
        if !(QoI_dict[QoI[i]][1] < x[i] < QoI_dict[QoI[i]][3])
            return false
        end
    end
    return true
end

function energy(x)
    # Call the Forward Map from the Package
    # It returns the 4 simulated trajectories
    tcan_p, tair_p, rh_p, co2_p = FM_SimLoop!(x, QoI, base_ctx, days, t_span, t_eval)
    
    # Calculate Log-Likelihood (Sum of Squares)
    logL = - sum((tcan_obs .- tcan_p).^2) / (error_std[1]^2) -
             sum((tair_obs .- tair_p).^2) / (error_std[2]^2) -
             sum((rh_obs   .- rh_p).^2)   / (error_std[3]^2) -
             sum((co2_obs  .- co2_p).^2)  / (error_std[4]^2)
    return -logL
end

# 8. --- Run MCMC ---
println("Starting MCMC for parameters: ", QoI)
n_params = length(QoI)
mcmc_obj = jtwalk(n=n_params, U=energy, Supp=PriorSupp)

# Initial points
x0  = [rand() * (QoI_dict[q][3] - QoI_dict[q][1]) + QoI_dict[q][1] for q in QoI]
xp0 = [rand() * (QoI_dict[q][3] - QoI_dict[q][1]) + QoI_dict[q][1] for q in QoI]

Run!(mcmc_obj, T=chain_length, x0=x0, xp0=xp0)

# 9. --- Plotting ---
function save_plots(mcmc, name, idx)
    p_chain = plot(mcmc.Output[:, idx], title="Chain: $name", label="")
    p_hist  = histogram(mcmc.Output[burn_in:end, idx], title="Posterior: $name", label="")
    vline!(p_hist, [QoI_dict[name][2]], color=:green, label="True/Prior Val")
    
    savefig(p_chain, "figures/$(name)_chain.png")
    savefig(p_hist, "figures/$(name)_posterior.png")
end

for (i, name) in enumerate(QoI)
    save_plots(mcmc_obj, name, i)
end