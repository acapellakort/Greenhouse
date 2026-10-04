# =============================================================================
# benchmark.jl  --  verify type stability, allocations, speed, and correctness
# =============================================================================
# Run from the GreenhouseSim folder:
#     julia --project=. test/benchmark.jl
#
# What to look for:
#   * @code_warntype: no red  ::Any / ::Union  on the hot variables.
#   * single rhs! call:  0 allocations, ~tens of ns.
#   * forward_map:        allocations in KB (not GB), fast.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using BenchmarkTools, Dates, InteractiveUtils

# --- paths to the proven 2.0 config + data (V2.2 sits beside invernaderos2.0)
INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))   # -> invernaderos root
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")
obs_path     = joinpath(INV, "invernaderos2.0", "observed_data.csv")

# --- build inputs ------------------------------------------------------------
start_date = DateTime(1998, 7, 11, 0, 0)
end_date   = DateTime(1998, 7, 13, 12, 0)     # cover the 2 simulated days

params   = load_params(climate_json, crop_json)
weather  = load_weather(meteo, start_date, end_date)
controls = load_measured_controls(obs_path)
base     = SimBase(params, weather, controls, datetime2unix(start_date))

days  = [0.0, 1.0]
tspan = (0.0, 86400.0)
teval = 0.0:3600.0:86400.0

# --- 1. single RHS call: type stability + allocations ------------------------
u   = copy(GreenhouseSim.DEFAULT_U0)
du  = similar(u)
ctx = GreenhouseSim.SimContext(params, weather, controls, datetime2unix(start_date), 0.0)

println("\n===== @code_warntype rhs!  (want: no red ::Any/::Union) =====")
@code_warntype rhs!(du, u, ctx, 0.0)

println("\n===== @btime single rhs! call  (target: 0 allocations) =====")
@btime rhs!($du, $u, $ctx, 0.0)

println("\n===== @btime single climate_rhs kernel =====")
@btime GreenhouseSim.climate_rhs(291.15, 296.15, 575.0, 1200.0,
    2.0, 300.0, 296.15, 272.75, 290.15, 273.15, 291.15, 3.0, 200.0, 834.7, 1000.0, 4.0,
    0,0,0,0,0,0,0,0,0,0,1, $params)

# --- 2. full forward map -----------------------------------------------------
QoI = ["beta2", "gamma3", "gamma4", "nu4"]
x   = [0.7, 275.0, 82.0, 0.01]

println("\n===== @btime forward_map (full 2-day solve, 4 params) =====")
@btime forward_map($x, $QoI, $base, $days, $tspan, $teval)

# --- 3. sanity / regression --------------------------------------------------
T1, T2, RH, CO2 = forward_map(x, QoI, base, days, tspan, teval)
println("\nlength(T2) = ", length(T2), "   (observed_data.csv has 51 rows)")
println("T2[1:3]  (K) = ", round.(T2[1:3], digits = 3))
println("RH[1:3]  (%) = ", round.(RH[1:3], digits = 3))
println("CO2[1:3]     = ", round.(CO2[1:3], digits = 3))
println("\nOK. Compare these against the 2.0 forward run for the regression check.")
