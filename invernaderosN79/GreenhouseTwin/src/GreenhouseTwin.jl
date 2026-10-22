"""
    GreenhouseTwin  (n79, flag-free)

The greenhouse climate + FvCB photosynthesis + cucumber crop model of the AiCU compartment (AGC-2018), in the
configuration selected as the baseline for the physics work: the V2.2 engine with the n79 flag set
(PRIORS=n79, SOLARC=1) resolved into one equation set, and Aaron Velez's generic canopy crop model (n105 set-up).
There are NO behaviour flags: every process below is always on. Documentation: docs/greenhouse_model_n79.tex.

Built from GreenhouseSim V2.2 at tag v2.2-flags-final (git commit 614a178); test/equivalence_n79.jl checks that
this code reproduces V2.2 with the n79 flags.

Climate state (per m2 floor):  u = [T1 canopy, T2 air, C1 CO2, V1 vapour pressure, e_int heating-error integral,
                                    T5 floor layer, Td deep soil, Tm internal mass]
Crop: daily step on measured climate (offline coupling): its LAI, respiration and fresh-mass series are climate inputs.
"""
module GreenhouseTwin

using OrdinaryDiffEq
using JSON
using CSV
using DataFrames
using Dates
using Statistics

include("parameters.jl")
include("inputs.jl")
include("psychrometrics.jl")
include("solar.jl")
include("photosynthesis.jl")
include("climate.jl")
include("simulate.jl")
include("agc_data.jl")
include("climate_calibration.jl")
include("crop.jl")
include("crop_measured.jl")

export load_params, update_params, param_meta, inferred_names, priors_of
export StepFun, ConstFun, WeatherInputs, ControlInputs
export Pws, rhf
export SolarFluxes, solar_column, column4, canopy_layer, layer_eff, roof_props, par_both_faces_factor
export assimilation, lai_eff
export SimContext, rhs!, model_terms, cover_temperature, NSTATE
export simdays, load_agc_season
export CLIMATE_QOI, CLIMATE_NOISE, climate_blocks, climate_energy, block_sse, season_eval, design_sets
export CucumberGeneric, cucumber_generic, CropState, Fruit, crop_step!, seasonal_sla, init_crop_state
export crop_run, crop_residuals, CROP_NOISE
export MeasuredSeason, load_measured_climate_par, load_fruit_keep, load_cropmanagement, load_production
export simulate_crop_measured, daily_assimilation_measured, crop_observations, crop_energy, CROP_QOI, crop_series

end # module
