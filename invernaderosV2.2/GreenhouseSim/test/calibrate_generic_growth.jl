# =============================================================================
# calibrate_generic_growth.jl -- crop calibration with Aaron's generic (Marcelis 1994) model (Note 69)
# =============================================================================
# Run from the GreenhouseSim folder:
#     MODE=centre julia --project=. test/calibrate_generic_growth.jl          (seconds: fit at the prior centre)
#     MODE=smoke  julia --project=. test/calibrate_generic_growth.jl          (~2000 steps: timing + sanity)
#     MODE=full SEED=1 julia --project=. test/calibrate_generic_growth.jl     (long run -- ask first)
#
# Driver: MEASURED AiCU air T and CO2; canopy PAR from the organisers' formula
#   PAR_in = PARout * transmission (0.65) * (1 - 0.98 BlackScr) * (1 - 0.25 EnScr) + 187 umol * lamps.
# Crop model: generic path (CUCUMBER_PARAMS = Marcelis 1994 via Aaron's port), with
#   node rate = flower appearance (per stem), M94 light factor on the daily PAR (gen_parlight = 1),
#   T clamped to the M94 domain (gen_clampT = 1), harvest when a fruit reaches 16 g DM (~400 g FW),
#   fruits set = nodes x stems x KEPT SHARE (grower policy, CropManagement 'Pruning') x carbon gate (abortion).
# Observations: weekly harvest (Production.csv), leaves per stem (N_leaves), fruit development time (FruitGrw).
#   The derived LAI is NOT used: it is built from assumed leaf area and stem density (Note 68).
# Inferred: gen_beta, gen_kappa, gen_A, gen_V, LAI_max, set_r_low, set_r_high, stem_density + 3 sigmas.
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, Dates, JTwalk, JLD2, Statistics, Printf, Random

const MODE  = get(ENV, "MODE", "centre"); MODE in ("centre", "smoke", "full") || error("MODE = centre | smoke | full")
const SEED  = parse(Int, get(ENV, "SEED", "1")); Random.seed!(SEED)
const NCHN  = parse(Int, get(ENV, "CHAIN", MODE == "smoke" ? "2000" : "150000"))
const BURN  = parse(Int, get(ENV, "BURN", "40000"))
const TRANS = parse(Float64, get(ENV, "TRANS", "0.65"))
const TEAM  = "AiCU"
const ZCLAMP = parse(Float64, get(ENV, "ZCLAMP", "0"))   # Note 70: 1 = fruits grow past stage 1, picked by size only

INV   = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA  = joinpath(INV, "invernaderos2.0", "data")
params = load_params(joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json"),
                     joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json"))
gp_base = load_growth_params(joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json"))
gp_base = update_params(gp_base, ["set_start_day", "harvest_W", "gen_clampT", "gen_parlight", "gen_beta", "gen_kappa",
                                  "gen_A", "gen_V", "gen_zclamp", "stem_density"], [15.0, 16.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, ZCLAMP, 2.6])

season = load_measured_climate_par(joinpath(DATA, TEAM, "Greenhouse_climate.csv"), joinpath(DATA, "meteo.csv");
                                   transmission = TRANS, lamp_umol = 187.0)
cm_csv = joinpath(DATA, TEAM, "CropManagement.csv")
keep   = load_fruit_keep(cm_csv; ndays = season.ndays)

# --- observations --------------------------------------------------------------------------------------------
days_obs, y_obs = load_production(joinpath(DATA, TEAM, "Production.csv"))
cum_at(dd, cum, d) = (d <= dd[1] ? float(cum[1]) : d >= dd[end] ? float(cum[end]) :
    let i = searchsortedlast(dd, d); dd[i] == d ? float(cum[i]) : cum[i] + (d - dd[i]) / (dd[i+1] - dd[i]) * (cum[i+1] - cum[i]) end)
weekly_incr(dd, cum, e) = diff([cum_at(dd, cum, x) for x in e])
const LAST = min(season.ndays - 1, days_obs[end])
const edges = let e = collect(0:7:LAST); e[end] < LAST && push!(e, LAST); e end
const Δy_obs = weekly_incr(days_obs, y_obs, edges)
cm_days, cm_N, _ = load_cropmanagement(cm_csv)
keepcm = cm_days .<= season.ndays - 1; cm_days = cm_days[keepcm]; cm_N = cm_N[keepcm]
cmdf = DataFrame(CSV.File(cm_csv))
fd_start = [Dates.value((Date(2018, 1, 1) + Day((Int(w) - 1) * 7)) - Date(2018, 8, 14)) for w in cmdf.weeks]
fd_val   = [x isa Number && !isnan(x) ? Float64(x) : NaN for x in cmdf.FruitGrw]
okfd     = .!isnan.(fd_val) .& (fd_start .>= 0) .& (fd_start .+ 6 .<= season.ndays - 1)
fd_start = fd_start[okfd]; fd_val = fd_val[okfd]
at_day(v, d) = v[min(d, length(v) - 1) + 1]
@printf("Observations: %d weekly yields, %d leaf counts, %d fruit-development weeks; kept share %.2f -> %.2f\n",
        length(Δy_obs), length(cm_N), length(fd_val), keep[1], keep[end])

function week_age(r, d0)     # n-weighted mean age of fruits harvested in days d0..d0+6
    n = 0.0; a = 0.0
    for d in d0:min(d0 + 6, length(r.day) - 1)
        hn = r.harv_n[d + 1]; hn > 0 && (n += hn; a += hn * r.harv_age[d + 1])
    end
    n > 0 ? a / n : NaN
end

# --- parameters ----------------------------------------------------------------------------------------------
QoI_G = ["gen_beta", "gen_kappa", "gen_A", "gen_V", "LAI_max", "set_r_low", "set_r_high"]   # stem_density fixed 2.6 (user, Note 71)
QoI_S = ["sigma_y", "sigma_leaf", "sigma_fd"]
QoI   = vcat(QoI_G, QoI_S); NG = length(QoI_G)
PR = Dict("gen_beta" => [0.5, 1.0, 1.5],  "gen_kappa" => [0.5, 1.0, 2.5], "gen_A" => [0.2, 1.0, 2.0],
          "gen_V" => [0.2, 1.0, 2.0],     "LAI_max" => [1.5, 2.5, 4.0],   "set_r_low" => [0.0, 0.3, 1.0],
          "set_r_high" => [0.2, 1.0, 2.5],
          "sigma_y" => [0.02, 0.3, 3.0],  "sigma_leaf" => [0.5, 3.0, 20.0], "sigma_fd" => [0.5, 3.0, 15.0])
centre = [PR[q][2] for q in QoI]

runmodel(z) = simulate_growth_measured(params, update_params(gp_base, QoI_G, z[1:NG]), season;
                                       LAI0 = 0.5, generic_params = CUCUMBER_PARAMS, fruit_keep = keep)
function resid(r)
    ry = Δy_obs .- weekly_incr(r.day, r.yield_FW, edges)
    rl = cm_N .- [at_day(r.node, d) for d in cm_days]
    rf = [isnan(week_age(r, d)) ? v : v - week_age(r, d) for (d, v) in zip(fd_start, fd_val)]   # no harvest = full miss
    ry, rl, rf
end
function energy(z)
    r = try runmodel(z) catch; return Inf end
    ry, rl, rf = resid(r); σy, σl, σf = z[NG + 1], z[NG + 2], z[NG + 3]
    length(ry) * log(σy) + sum(abs2, ry) / (2σy^2) + length(rl) * log(σl) + sum(abs2, rl) / (2σl^2) +
    length(rf) * log(σf) + sum(abs2, rf) / (2σf^2)
end
supp(z) = all(PR[q][1] < z[k] < PR[q][3] for (k, q) in enumerate(QoI)) && z[7] > z[6] + 0.05

function report(z, tag)
    r = runmodel(z); ry, rl, rf = resid(r)
    @printf("\n  %s: final yield %.2f kg FW/m2 (obs %.2f) | weekly-yield RMSE %.3f | leaves RMSE %.2f (final %.1f vs %.1f) | fruit-age RMSE %.2f d\n",
            tag, r.yield_FW[end], cum_at(days_obs, y_obs, LAST), sqrt(mean(abs2, ry)), sqrt(mean(abs2, rl)),
            at_day(r.node, cm_days[end]), cm_N[end], sqrt(mean(abs2, rf)))
    @printf("  fruits harvested %.1f /m2, mean fruit FW %.0f g; final LAI %.2f; mean Pg %.1f g CH2O/m2/d (d0-48 %.1f, d49-114 %.1f)\n",
            sum(r.harv_n), sum(r.harv_n) > 0 ? 1000 * r.yield_FW[end] / sum(r.harv_n) : NaN, r.LAI[end], mean(r.Pg),
            mean(r.Pg[1:49]), mean(r.Pg[50:end]))
    println("  week start | fruit dev. time obs / model (d)")
    for (d, v) in zip(fd_start, fd_val); @printf("  %10d | %5.1f / %5.1f\n", d, v, week_age(r, d)); end
end

println("\nGeneric (Aaron/M94) crop calibration -- mode=$MODE seed=$SEED transmission=$TRANS zclamp=$ZCLAMP")
report(centre, "prior centre")
@printf("  energy at centre %.1f\n", energy(centre))
MODE == "centre" && exit()

obj = jtwalk(n = length(QoI), U = energy, Supp = supp)
jit() = (for _ in 1:1000; x = centre .* (1 .+ 0.05 .* randn(length(QoI))); supp(x) && return x; end; centre)
t0 = time(); Run!(obj, T = NCHN, x0 = jit(), xp0 = jit())
@printf("\n%d steps in %.0f s (%.2f s/step)\n", NCHN, time() - t0, (time() - t0) / NCHN)
postdir = joinpath(@__DIR__, "posteriors"); mkpath(postdir)
fn = joinpath(postdir, "$(TEAM)_growth_generic$(ZCLAMP > 0.5 ? "_zc" : "")_seed$(SEED)$(MODE == "smoke" ? "_smoke" : "").jld2")
JLD2.@save fn chain=obj.Output QoI=QoI burn_in=BURN transmission=TRANS
println("saved ", fn)
b = MODE == "smoke" ? max(1, NCHN ÷ 2) : BURN
post = obj.Output[b:end, 1:end-1]; pm = vec(mean(post, dims = 1)); ps = vec(std(post, dims = 1))
println("\nposterior (mean ± sd):"); for (k, q) in enumerate(QoI); @printf("  %-13s %8.4f ± %.4f\n", q, pm[k], ps[k]); end
report(pm, "posterior mean")
