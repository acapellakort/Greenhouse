"""
    GreenhouseControl  (V3.0)

Control / agent layer for the greenhouse digital twin (project goal #2).

REUSES the frozen V2.2 model engine (`GreenhouseSim`: climate physics + FvCB
photosynthesis + cucumber growth) by `include` — it does NOT copy that code, so
there is a single source of truth for the model. This package adds only the new
control layer (and, later, the RL environment and agent):

  * a "climate computer": setpoint schedules + PI/proportional controllers that
    turn setpoints into the actuators U1..U12 (control on air temp T2);
  * a closed-loop RHS (`rhs_control!`) driven by those setpoints;
  * (Stage 2) coupling to the growth model day-by-day;
  * (Stage 3) an RL environment: action = setpoints, reward = net profit.
"""
module GreenhouseControl

using OrdinaryDiffEq
using DataInterpolations
using Dates
using Random

# --- frozen V2.2 engine (include, do NOT copy) -------------------------------
include(joinpath(@__DIR__, "..", "..", "..",
                 "invernaderosV2.2", "GreenhouseSim", "src", "GreenhouseSim.jl"))
using .GreenhouseSim

include("climate_computer.jl")
include("twin.jl")
include("economics.jl")
include("weather_episodes.jl")
include("rl_env.jl")

# control layer
export Setpoints, daynight_setpoints, PIGains, vent_control, screen_control, humidity_vent
export ControlContext, rhs_control!, run_control
# coupled twin (Stage 2)
export simulate_twin, day_resources, density_schedule
export season_economics, print_economics, PRICE_CUKE
export TwinEnv, env_reset!, env_step!, observe, decode_setpoints, ACT_LO, ACT_HI, ACT_NAMES
export WeatherBank, valid_start_years, sample_window, SeasonForecast, forecast_row
# convenience re-exports from the engine
export load_params, load_weather, WeatherInputs
export simulate_growth_measured, load_growth_params, update_params
export init_growth_state, grow!, GrowthState

end # module
