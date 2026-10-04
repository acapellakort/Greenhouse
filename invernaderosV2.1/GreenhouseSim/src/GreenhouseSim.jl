module GreenhouseSim

using Reexport
using DifferentialEquations
@reexport using DataInterpolations
using Unitful
using JSON
using Parameters
using Dates
using DataFrames
using CSV

# --- Include Files (ONLY ONCE EACH) ---
include("utils/Parameters.jl")
include("utils/Controls.jl")
include("models/climate_functionsV3_NoUnits.jl")
include("models/cropphoto_functionsV3_NoUnits.jl")
include("models/temp_controlsV3.jl")
include("models/physics_ode.jl")

# --- Export everything needed by run_inference.jl ---
export FM_SimLoop!, rhsPLL_PID_sim!
export ModelParams, load_model_params, merge_params, unit_mapping
export ControlSet, MeasuredControls, load_controls, load_control_interpolations, get_val
export SimulationContext

# --- SimulationContext Definition ---
@with_kw struct SimulationContext
    params::ModelParams
    controls_json::ControlSet
    controls_csv::MeasuredControls
    weather::Dict{Symbol, Any}
    start_time_unix::Float64
    start_date::DateTime
end
# --- 5. Forward Map / Inference Logic ---
"""
    FM_SimLoop!(x, theta_names, base_ctx, days, t_span, t_eval)
This is the MCMC loop. It creates a new context for every 'x' 
proposed by the sampler without using global variables.
"""
function FM_SimLoop!(x, theta_names, base_ctx::SimulationContext, days, t_span, t_eval)
    # Update parameters locally for this MCMC step
    new_vals = merge(base_ctx.params.vals, NamedTuple{Tuple(Symbol.(theta_names))}(x))
    new_params = ModelParams(new_vals, base_ctx.params.units, base_ctx.params.meta)
    
    # Reconstruct context with updated parameters
    ctx = SimulationContext(
        params = new_params,
        controls_json = base_ctx.controls_json,
        controls_csv = base_ctx.controls_csv,
        weather = base_ctx.weather,
        start_time_unix = base_ctx.start_time_unix,
        start_date = base_ctx.start_date
    )

    # Initial conditions (can be passed in or set here)
    u0 = [18+273.15, 23+273.15, 575, 1200, 0.0] 

    # Run the days
    for day_i in days
        # We pass the ctx (which contains all params/weather/controls) as 'p'
        prob = ODEProblem(rhsPLL_PID_sim!, u0, t_span, ctx)
        sol = solve(prob, Tsit5(), saveat=t_eval)
        u0 = sol[:, end] # Carry over state to next day
    end
    
    # Return simulated outputs...
    # (Add return logic based on your Likelihood requirements)
end

end # module GreenhouseSim