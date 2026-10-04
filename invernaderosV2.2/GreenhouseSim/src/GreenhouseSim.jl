"""
    GreenhouseSim  (V2.2)

Type-stable, allocation-free rewrite of the greenhouse climate + FvCB
photosynthesis forward map, built from the proven 2.0 physics.

Design goals
------------
1. SPEED: no non-`const` globals anywhere on the hot path. Every physical
   constant is carried in a single concretely-typed parameter `NamedTuple`
   (`Float64` values) and read with literal destructuring `(; a, b) = p`,
   which the compiler fully specializes. Weather and controls are concrete
   `NamedTuple`s of callable interpolants (never `Dict{Symbol,Any}`).

2. INFERENCE FLEXIBILITY (the key 2.0 feature, kept without globals):
   choose *any* subset of parameters to infer, by name, at run time.
   `update_params(base, names, x)` returns a new concrete parameter
   `NamedTuple` with exactly those names overridden — the same freedom the
   old globals gave, but the (tiny) dynamic cost happens once per MCMC
   proposal, OUTSIDE the ODE inner loop (a function barrier).

3. GENERIC DRY-MATTER PARTITIONING (default-off):
   Pass `generic_params = CUCUMBER_PARAMS` or `TOMATO_PARAMS` to `grow!`
   / `grow_cohort!` / `simulate_growth` / `simulate_cohort_measured` to
   activate the Marcelis (1994) / Heuvelink (1996) Bell-curve sink and
   affine appearance rate. Default `nothing` preserves every V2.2 result.

Public API
----------
    load_params(climate_json, crop_json) -> NamedTuple      # all constants
    update_params(base, names, x)         -> NamedTuple      # override by name
    load_weather(csv, start, stop)        -> WeatherInputs
    load_measured_controls(csv)           -> ControlInputs
    SimContext(params, weather, controls, weather_t0, controls_t0)
    rhs!(du, u, ctx, t)                                     # ODE right-hand side
    forward_map(x, qoi, base_ctx, days, tspan, teval; u0)  # -> (T1, T2, RH, CO2)

See test/benchmark.jl and test/run_inference.jl for usage.
"""
module GreenhouseSim

using OrdinaryDiffEq
using DataInterpolations
using JSON
using Dates
using DataFrames
using CSV

include("parameters.jl")
include("weather.jl")
include("controls.jl")
include("climate.jl")
include("photosynthesis.jl")
include("generic_growth.jl")   # GenericParameters, CUCUMBER_PARAMS, TOMATO_PARAMS
include("forward_map.jl")
include("growth.jl")
include("measured.jl")
include("growth_cohort.jl")

# Parameters
export load_params, update_params
# Inputs
export WeatherInputs, load_weather, prepare_meteo, weather_from_df
export ControlInputs, load_measured_controls
# Simulation
export SimBase, SimContext, rhs!, simulate, forward_map
# Model kernels (exported for benchmarking / testing)
export climate_rhs, assimilation, Pws, rhf, VPDf, floor_balance
# Generic dry-matter partitioning (Marcelis 1994 / Heuvelink 1996)
export GenericParameters, CUCUMBER_PARAMS, TOMATO_PARAMS
export generic_appearance, generic_development_rate, generic_sink, generic_vegetation
export allocate_daily
# M94 / H96 reference functions (for validation)
export cucumber_appearance_m94, cucumber_sink_m94, cucumber_vegetation_m94
export tomato_appearance_h96, tomato_development_rate_h96, tomato_sink_h96, tomato_vegetation_h96
# Growth model (cucumber source-sink)
export load_growth_params, Fruit, GrowthState, init_growth_state
export rhs_dyn!, day_assimilation, grow!, simulate_growth
export LeafCohort, CohortState, init_cohort_state, grow_cohort!, simulate_cohort_measured
export daily_assimilation_layered, canopy_light_profile
# Measured-climate calibration support
export MeasuredSeason, load_measured_climate, simulate_growth_measured
export daily_assimilation_measured, load_production, load_cropmanagement

end # module
