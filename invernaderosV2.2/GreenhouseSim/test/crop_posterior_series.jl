# =============================================================================
# crop_posterior_series.jl -- daily crop series from the n78 crop posterior, as inputs to the n79 climate run (Note 79)
# =============================================================================
#     julia --project=. test/crop_posterior_series.jl            (~10 s; writes test/posteriors/AiCU_crop_series_n78.csv)
#
# Runs the calibrated Aaron-V2 crop (same set-up as calibrate_generic_growth.jl: measured T/CO2, PAR at 0.65, kept-fruit
# policy, stem density 2.6, picking at 400 g FW) at the POOLED posterior mean of the n78 seeds (burn-in 40000), and writes
#   LAI_floor   = LAI x 0.8  (crop LAI per growing area -> per floor area; growing area 0.8 of the floor)
#   CropResp20  = non-leaf maintenance + growth respiration at 20 C, mg CO2 m-2 floor s-1:
#                 0.8 * [ 2^((20-25)/10) (k_m_stem W_stem + k_m_root W_root + k_m_fruit W_fruit) + (asrq-1) max(dDM,0) ]
#                 * 44/30 * 1000/86400        (leaf maintenance is excluded: it is the FvCB dark respiration Rd)
# The climate model applies Q10 = 2 on the canopy temperature (forward_map.jl, flag co2_resp).
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Dates, JLD2, Statistics, Printf

const TRANS = 0.65; const TEAM = "AiCU"; const BURN = 40000
const CROPTAG = get(ENV, "CROPTAG", "n78")       # Note 94: n94 = canopy model, everything already per m2 floor
const CANOPY  = CROPTAG in ("n94", "n105")      # Note 105: n105 = n94 + par_below
const FLOOR   = CANOPY ? 1.0 : 0.8
INV  = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA = joinpath(INV, "invernaderos2.0", "data")
params = load_params(joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json"),
                     joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json"))
CROPTAG == "n105" && (params = update_params(params, ["par_below", "tau1"], [1.0, 0.65 / 0.9]))
gp_base = load_growth_params(joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json"))
gp_base = update_params(gp_base, ["set_start_day", "harvest_W", "gen_clampT", "gen_parlight", "gen_beta", "gen_kappa",
                                  "gen_A", "gen_V", "gen_zclamp", "stem_density"], [15.0, 16.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, CANOPY ? 3.6 : 2.6])
CANOPY && (gp_base = update_params(gp_base, ["gen_canopy", "gen_doy0", "gen_fL", "gen_flower_node"], [1.0, 226.0, 0.59, 5.0]))
season = load_measured_climate_par(joinpath(DATA, TEAM, "Greenhouse_climate.csv"), joinpath(DATA, "meteo.csv");
                                   transmission = TRANS, lamp_umol = 187.0)
keep = load_fruit_keep(joinpath(DATA, TEAM, "CropManagement.csv"); ndays = season.ndays, as_removal = CANOPY)   # Note 97

# pooled posterior mean of the n78 seeds
chains = Matrix{Float64}[]; QoI = String[]
for s in 1:2
    fn = joinpath(@__DIR__, "posteriors", CANOPY ? "$(TEAM)_growth_canopy_$(CROPTAG)_seed$(s).jld2" : "$(TEAM)_growth_generic_n78_seed$(s).jld2")
    isfile(fn) || (println("missing ", fn); continue)
    d = JLD2.load(fn); push!(chains, d["chain"][BURN:end, 1:end-1]); global QoI = d["QoI"]
end
isempty(chains) && error("no n78 crop chains found")
pm = vec(mean(vcat(chains...), dims = 1)); ix(q) = findfirst(==(q), QoI)
println("n78 pooled posterior mean ($(length(chains)) seeds):")
for (k, q) in enumerate(QoI); @printf("  %-11s %.4g\n", q, pm[k]); end

const V2 = CUCUMBER_PARAMS_V2
function genpar(W_pot, E, beta1, V20)
    AT = W_pot * V2.s * (V2.p - 1) / E
    GenericParameters(s = V2.s, p = V2.p, peak = V2.peak, A0 = -10.0 * AT, AT = AT, kappa = 1.0 / E, alpha = 1.0,
                      base = 10.0, beta0 = -beta1 * (0.75 / 0.09), beta1 = beta1, V20 = V20, VT = V2.VT,
                      clock_domain = V2.clock_domain, appearance_domain = V2.appearance_domain,
                      sink_domain = V2.sink_domain, vegetation_domain = V2.vegetation_domain)
end
QoI_GP = CANOPY ? ["k_m_leaf", "k_m_stem", "k_m_fruit", "gen_leaf_target", "gen_N0", "set_r_low", "set_r_high", "DMC_fruit"] :
                  ["k_m_leaf", "k_m_stem", "k_m_fruit", "LAI_max", "set_r_low", "set_r_high", "DMC_fruit"]
gp = update_params(gp_base, vcat(QoI_GP, ["harvest_W"]), vcat([pm[ix(q)] for q in QoI_GP], (CANOPY ? 323.0 : 400.0) * pm[ix("DMC_fruit")]))
r = simulate_growth_measured(update_params(params, ["J_max"], [pm[ix("J_max")]]), gp, season; LAI0 = CANOPY ? 0.25 : 0.5,
                             fruit_keep = keep, generic_params = genpar(pm[ix("W_pot")], pm[ix("E_fruit")], pm[ix("beta1")], pm[ix("V20")]))

DM   = r.W_leaf .+ r.W_stem .+ r.W_root .+ r.W_fruit .+ r.yield_DW
dDM  = max.(vcat(0.0, diff(DM)), 0.0)
Rm20 = 2.0^((20 - 25) / 10) .* (gp.k_m_stem .* r.W_stem .+ gp.k_m_root .* r.W_root .+ gp.k_m_fruit .* r.W_fruit)
Rg   = (gp.asrq - 1) .* dDM
conv = 44 / 30 * 1000 / 86400
CR20 = FLOOR .* (Rm20 .+ Rg) .* conv
out  = DataFrame(day = r.day, LAI_floor = FLOOR .* r.LAI, CropResp20 = CR20, Rm20_gCH2O = Rm20, Rg_gCH2O = Rg,
                 Pg_gCH2O = r.Pg, yield_FW = r.yield_FW, node = r.node, W_leaf = r.W_leaf, W_stem = r.W_stem,
                 W_root = r.W_root, W_fruit = r.W_fruit, yield_DW = r.yield_DW, harv_n = r.harv_n, LAI = r.LAI,
                 # Note 91: standing fresh mass per m2 FLOOR of stems + roots + fruits on the plant (leaves are in the canopy
                 # heat capacity alpha1*LAI). Dry-matter contents: stem 0.06, root 0.07 (assumed), fruit = inferred DMC_fruit.
                 FW_mass_kg = FLOOR .* (r.W_stem ./ 0.06 .+ r.W_root ./ 0.07 .+ r.W_fruit ./ pm[ix("DMC_fruit")]) ./ 1000)
fo = joinpath(@__DIR__, "posteriors", "$(TEAM)_crop_series_$(CROPTAG).csv"); CSV.write(fo, out)
@printf("\nfinal yield %.2f kg FW/m2; LAI_floor mean %.2f (max %.2f)\n", r.yield_FW[end], mean(out.LAI_floor), maximum(out.LAI_floor))
@printf("CropResp20 mean %.4f mg CO2/m2/s (maint %.4f, growth %.4f); = %.1f%% of mean Pg in CO2 terms\n",
        mean(CR20), mean(FLOOR .* Rm20 .* conv), mean(FLOOR .* Rg .* conv), 100 * mean(Rm20 .+ Rg) / mean(r.Pg))
println("weekly (day, LAI_floor, CropResp20):")
for d in 0:14:length(r.day)-1; @printf("  %4d  %5.2f  %.4f\n", r.day[d+1], out.LAI_floor[d+1], CR20[d+1]); end
@printf("standing fresh mass (stems+roots+fruits, floor basis): day 0 %.2f, day 30 %.2f, end %.2f kg/m2\n",
        out.FW_mass_kg[1], out.FW_mass_kg[min(31, end)], out.FW_mass_kg[end])
println("wrote ", fo)
