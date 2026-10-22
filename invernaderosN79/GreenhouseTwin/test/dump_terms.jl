# =============================================================================
# dump_terms.jl -- every term of model_terms at random states and times (for tools/check_tex.py)
# =============================================================================
#     julia --project=. test/dump_terms.jl            -> test/terms_dump.csv  (synthetic inputs, no data files)
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin, Random, CSV, DataFrames
const GT = GreenhouseTwin
p = load_params(joinpath(@__DIR__, "..", "configfiles", "climate_n79.json"))
T0 = 1.5342e9; n = 2 * 288 + 1; ts = [T0 + 300.0 * (k - 1) for k in 1:n]; tr = ts .- T0
rng = MersenneTwister(11)
rv(lo, hi) = lo .+ (hi - lo) .* rand(rng, n)
inp = Dict(:Iglob => rv(0, 900), :Tsky => rv(255, 290), :Tout => rv(275, 300), :Wind => rv(0, 8), :VPout => rv(500, 1800),
           :LAI => rv(0.3, 4), :EnScr => rv(0, 1), :BlackScr => rv(0, 1), :FanPad => rv(0, 1), :MechCool => rv(0, 1),
           :Heater => rv(0, 1), :Shade => rv(0, 1), :SideVent => rv(0, 1), :ForcedVent => rv(0, 1), :VentLee => rv(0, 1),
           :VentWind => rv(0, 1), :Fog => rv(0, 1), :CO2Valve => rv(0, 1), :HeatSP => rv(288, 294), :Lamps => rv(0, 1),
           :Tpipe => rv(285, 330), :CropResp => rv(0, 0.01), :CropFW => rv(1, 8))
S(k) = StepFun(k in (:Iglob, :Tsky, :Tout, :Wind, :VPout, :LAI) ? ts : tr, inp[k])
w = WeatherInputs(S(:Iglob), S(:Tsky), S(:Tout), S(:Wind), S(:VPout), ConstFun(834.7), ConstFun(276.0), S(:LAI))
c = ControlInputs(S(:EnScr), S(:BlackScr), S(:FanPad), S(:MechCool), S(:Heater), S(:Shade), S(:SideVent), S(:ForcedVent),
                  S(:VentLee), S(:VentWind), S(:Fog), S(:CO2Valve), S(:HeatSP), S(:Lamps), S(:Tpipe), S(:CropResp), S(:CropFW))
ctx = SimContext(p, w, c, T0, 0.0)
rows = Vector{Dict{String,Float64}}()
for k in 1:400
    t = 300.0 * (rand(rng, 0:575)) + 150.0                    # mid-interval: inputs unambiguous
    T2 = 283 + 22rand(rng); RH = 30 + 72rand(rng)            # includes supersaturated air
    u = [T2 + 4randn(rng), T2, 350 + 1300rand(rng), RH / 100 * Pws(T2), 4000rand(rng) * (rand(rng) < 0.3 ? -1 : 1),
         T2 - 6 + 9rand(rng), 283 + 12rand(rng), T2 - 4 + 8rand(rng)]
    m = model_terms(u, ctx, t)
    d = Dict{String,Float64}("t" => t)
    for (j, nm) in enumerate(("T1", "T2", "C1", "V1", "eint", "T5", "Td", "Tm")); d[nm] = u[j]; end
    for (key, f) in inp; d["in_" * String(key)] = f[clamp(searchsortedlast(tr, t), 1, n)]; end
    for (key, val) in pairs(m)
        if val isa Real
            d[String(key)] = Float64(val)
        elseif val isa GT.SolarFluxes
            for fn in fieldnames(GT.SolarFluxes); d["so_" * String(fn)] = getfield(val, fn); end
        end
    end
    push!(rows, d)
end
cols = sort(collect(keys(rows[1])))
df = DataFrame([c => [r[c] for r in rows] for c in cols])
CSV.write(joinpath(@__DIR__, "terms_dump.csv"), df)
println("wrote test/terms_dump.csv: ", size(df))
