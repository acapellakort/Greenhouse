# =============================================================================
# diagnose_crop_bias_sensitivity.jl -- how much climate bias can the crop model tolerate? (Note 45b)
# =============================================================================
# Run from the GreenhouseSim folder (growth-model forward runs only, ~1 min, no climate model, no MCMC):
#     julia --project=. test/diagnose_crop_bias_sensitivity.jl
#
# Crop model: Note 21 posterior mean (calibrate_growth.jl), driven by MEASURED CO2 and light and by the
# measured air T plus a controlled bias. Day = canopy light sample > 0, night otherwise (5-min samples).
# Scenarios
#   uniform   bias added to every sample: ±0.25, ±0.5, ±1, ±2 K
#   night     bias only at night: ±0.5, ±1 K
#   day       bias only by day:   ±0.5, ±1 K
#   pattern   the actual per-day night/day bias of a climate model (from AiCU_season_daily.csv, if present):
#             N33 open loop, N38 closed loop, N39 what-if
# Outputs vs the unbiased run: final node number, final yield, yield on days 42-114 (holdout window of
# Note 31), final LAI, weekly-yield RMSE against observations. Then the bias that keeps |yield error| < 2 % / 5 %.
# Note: this crop model does not use humidity, so the VP criterion cannot be derived here.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = "115"
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))
using DelimitedFiles

gp21 = update_params(gp_base, QoI_G, centre[IG])
println("\nCrop parameters (Note 21 posterior mean): ", join(["$q=$(round(v, sigdigits=3))" for (q, v) in zip(QoI_G, centre[IG])], "  "))

# day/night flag per growth sample (light row of the sample matrix)
isday = [samp_meas[d][3, :] .> 0.0 for d in 1:N_DAYS]

function run_bias(fday, fnight)            # fday(d), fnight(d): bias in K on day d (0-based)
    for d in 1:N_DAYS
        M = work[d]; M0 = samp_meas[d]
        @inbounds for k in axes(M, 2)
            M[1, k] = M0[1, k] + (isday[d][k] ? fday(d - 1) : fnight(d - 1))
        end
    end
    simulate_growth_measured(params, gp21, season_B; LAI0 = 0.5)
end
yield_at(r, d) = r.yield_FW[min(d, length(r.yield_FW) - 1) + 1]
function metrics(r)
    wk = weekly_incr(r.day, r.yield_FW, edges)
    (node = r.node[end], yld = r.yield_FW[end], yhold = yield_at(r, 114) - yield_at(r, 42), lai = r.LAI[end],
     wrmse = sqrt(mean(abs2, wk .- Δy_obs)))
end

ref0 = metrics(run_bias(d -> 0.0, d -> 0.0))
@printf("Unbiased (measured T): nodes %.1f, yield %.2f kg/m2, yield d42-114 %.2f, LAI %.2f, weekly-yield RMSE vs obs %.3f\n",
        ref0.node, ref0.yld, ref0.yhold, ref0.lai, ref0.wrmse)
@printf("Observed: final yield %.2f kg/m2 (from Production.csv)\n", y_obs[end])

scen = Any[]
for b in (-2.0, -1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 2.0); push!(scen, ("uniform $(b > 0 ? "+" : "")$(b) K", b, b, d -> b, d -> b)); end
for b in (-1.0, -0.5, 0.5, 1.0); push!(scen, ("night only $(b > 0 ? "+" : "")$(b) K", 0.0, b, d -> 0.0, d -> b)); end
for b in (-1.0, -0.5, 0.5, 1.0); push!(scen, ("day only $(b > 0 ? "+" : "")$(b) K", b, 0.0, d -> b, d -> 0.0)); end

csvf = joinpath(FIGDIR, "AiCU_season_daily.csv")
if isfile(csvf)
    raw, hdr = readdlm(csvf, ',', header = true); hdr = vec(String.(hdr))
    colv(nm) = (j = findfirst(==(nm), hdr); j === nothing ? nothing : [x isa Number ? Float64(x) : NaN for x in raw[:, j]])
    for (k, lab) in (("OL33", "pattern N33 open loop"), ("CL38", "pattern N38 closed loop"), ("NV39", "pattern N39 what-if"))
        nb = colv("$(k)_night_bias"); db = colv("$(k)_day_bias")
        (nb === nothing || db === nothing) && continue
        nb = replace(nb, NaN => 0.0); db = replace(db, NaN => 0.0)
        push!(scen, (lab, mean(db), mean(nb), d -> db[min(d, length(db) - 1) + 1], d -> nb[min(d, length(nb) - 1) + 1]))
    end
else
    println("(no $csvf -- run test/full_season_report.jl first to add the model bias patterns)")
end

rows = Any[["scenario", "day_bias_K", "night_bias_K", "nodes", "dnodes_pct", "yield", "dyield_pct", "yield_42_114",
            "dyhold_pct", "LAI", "weekly_yield_RMSE"]]
println("\n", "="^128)
@printf("  %-26s | %6s %6s | %7s %7s | %7s %7s | %7s %7s | %5s | %7s\n", "scenario", "day", "night", "nodes", "Δ%",
        "yield", "Δ%", "y42-114", "Δ%", "LAI", "wkRMSE")
println("  ", "-"^126)
@printf("  %-26s | %6s %6s | %7.1f %7s | %7.2f %7s | %7.2f %7s | %5.2f | %7.3f\n", "unbiased (measured T)", "", "",
        ref0.node, "", ref0.yld, "", ref0.yhold, "", ref0.lai, ref0.wrmse)
pct(a, b) = 100 * (a - b) / b
for (nm_, db, nb, fd, fn) in scen
    m = metrics(run_bias(fd, fn))
    @printf("  %-26s | %+6.2f %+6.2f | %7.1f %+7.1f | %7.2f %+7.1f | %7.2f %+7.1f | %5.2f | %7.3f\n", nm_, db, nb,
            m.node, pct(m.node, ref0.node), m.yld, pct(m.yld, ref0.yld), m.yhold, pct(m.yhold, ref0.yhold), m.lai, m.wrmse)
    push!(rows, Any[nm_, db, nb, m.node, pct(m.node, ref0.node), m.yld, pct(m.yld, ref0.yld), m.yhold,
                    pct(m.yhold, ref0.yhold), m.lai, m.wrmse])
end
writedlm(joinpath(FIGDIR, "AiCU_crop_bias_sensitivity.csv"), permutedims(hcat(rows...)), ',')

# tolerance from the uniform scenarios: linear sensitivity of yield and nodes per K, then |bias| for 2 % / 5 %
Urows = [r for r in rows[2:end] if startswith(r[1], "uniform")]
bs = Float64[r[2] for r in Urows]; dy = Float64[r[7] for r in Urows]; dn = Float64[r[5] for r in Urows]
sy = sum(bs .* dy) / sum(bs .^ 2); sn = sum(bs .* dn) / sum(bs .^ 2)
@printf("\n  Linear sensitivity (uniform bias): yield %+.2f %% per K, nodes %+.2f %% per K\n", sy, sn)
@printf("  Bias that keeps |yield error| below 2 %%: %.2f K;  below 5 %%: %.2f K\n", 2 / abs(sy), 5 / abs(sy))
@printf("  Bias that keeps |node error|  below 2 %%: %.2f K;  below 5 %%: %.2f K\n", 2 / abs(sn), 5 / abs(sn))
println("""
  How to read it:
   - The tolerance lines turn a crop-forecast requirement (e.g. yield within 2 % or 5 %) into the climate bias
     acceptance criterion for Note 45 (replacing the ±0.5 K guess).
   - Compare night-only and day-only rows: if the crop model is less sensitive at night, the night criterion can
     be looser than the day criterion, or vice versa.
   - The 'pattern' rows show what the actual climate-model errors of Notes 33/38/39 would do to the forecast.
  Table written to test/figures/AiCU_crop_bias_sensitivity.csv""")
