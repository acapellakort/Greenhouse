# =============================================================================
# benchmark.jl -- type stability, allocations and speed of rhs! and of one simulated day (AiCU inputs)
# =============================================================================
#     CROPTAG=n94 julia --project=. test/benchmark.jl
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin
using BenchmarkTools, InteractiveUtils
const GT = GreenhouseTwin
INV = normpath(joinpath(@__DIR__, "..", "..", ".."))
tag = get(ENV, "CROPTAG", "n94")
csn = joinpath(@__DIR__, "posteriors", "AiCU_crop_series_$(tag).csv")
cs  = isfile(csn) ? csn : joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "test", "posteriors", "AiCU_crop_series_$(tag).csv")
p = load_params(joinpath(@__DIR__, "..", "configfiles", "climate_n79.json"))
S = load_agc_season(INV; crop_series = cs)
u = copy(S.u0); du = similar(u)
ctx = SimContext(p, S.weather, S.controls, S.t0 + 86400.0 * 50, 86400.0 * 50)
println("\n===== @code_warntype rhs! (want no ::Any / ::Union) =====")
@code_warntype rhs!(du, u, ctx, 43200.0)
println("\n===== rhs! (target: 0 allocations) =====")
@btime rhs!($du, $u, $ctx, 43200.0)
println("\n===== one day (reltol 1e-4) =====")
@btime simdays($p, $S, [50.0], $u, 1e-4)
