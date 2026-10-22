# =============================================================================
# crop_series.jl -- daily crop series (climate inputs) from a crop posterior
# =============================================================================
# Run from invernaderosN79/GreenhouseTwin:
#   SEEDS=2 BURN=150000 julia --project=. test/crop_series.jl
# Reads crop chains (QoI names as in crop_n105.json) from test/posteriors/AiCU_crop_n105_twin_seed<s>.jld2, or, with
# SOURCE=v22, the V2.2 chains AiCU_growth_canopy_n105_seed<s>.jld2 (QoI beta1 -> app_beta1). Runs the crop at the
# pooled posterior mean and writes test/posteriors/AiCU_crop_series_n105.csv (LAI_floor, CropResp20, FW_mass_kg, ...).
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin
using CSV, DataFrames, JLD2, Printf, Statistics
const GT = GreenhouseTwin

const SEEDS  = parse.(Int, split(get(ENV, "SEEDS", "2"), ","))
const BURN   = parse(Int, get(ENV, "BURN", "150000"))
const SOURCE = get(ENV, "SOURCE", "v22")          # v22 = V2.2 n105 chains (until the GreenhouseTwin chains exist)
const INV    = normpath(joinpath(@__DIR__, "..", "..", ".."))
const POST   = joinpath(@__DIR__, "posteriors"); mkpath(POST)
const pc  = load_params(joinpath(@__DIR__, "..", "configfiles", "climate_n79.json"))
const gp0 = load_params(joinpath(@__DIR__, "..", "configfiles", "crop_n105.json"))
const O   = crop_observations(INV; transmission = gp0.transmission, lamp_umol = gp0.lamp_umol)

chains = Matrix{Float64}[]; QoI = String[]
for s in SEEDS
    fn = SOURCE == "v22" ? joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "test", "posteriors", "AiCU_growth_canopy_n105_seed$(s).jld2") :
                           joinpath(POST, "AiCU_crop_n105_twin_seed$(s).jld2")
    isfile(fn) || (println("missing ", fn); continue)
    d = JLD2.load(fn); push!(chains, d["chain"][BURN:end, 1:end-1]); global QoI = replace(String.(d["QoI"]), "beta1" => "app_beta1")
end
isempty(chains) && error("no crop chains found")
pm = vec(mean(vcat(chains...), dims = 1)); ix(q) = findfirst(==(q), QoI)
println("posterior mean, seeds $(SEEDS), burn-in $(BURN):"); for (k, q) in enumerate(QoI); @printf("  %-16s %.4g\n", q, pm[k]); end
z = vcat([pm[ix(q)] for q in GT.CROP_QOI], [pm[ix(q)] for q in GT.CROP_NOISE])
r  = GT.crop_run(z, pc, gp0, O)
gp = update_params(gp0, GT.CROP_QOI[1:end-1], z[1:length(GT.CROP_QOI)-1])
out = crop_series(r, gp)
fo = joinpath(POST, "AiCU_crop_series_n105.csv"); CSV.write(fo, out)
@printf("\nfinal yield %.2f kg FW m-2; LAI mean %.2f (max %.2f); CropResp20 mean %.4f mg CO2 m-2 s-1\n",
        r.yield_FW[end], mean(out.LAI_floor), maximum(out.LAI_floor), mean(out.CropResp20))
println("wrote ", fo)
