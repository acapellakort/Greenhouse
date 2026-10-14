# =============================================================================
# diagnose_crop_n54.jl -- crop model driven by the n54 climate model instead of measured air T (Note 57; forward)
# =============================================================================
# Run from the GreenhouseSim folder (~2 min, no MCMC):
#     PRIORS=n54 julia --project=. test/diagnose_crop_n54.jl
#
# Crop model: Note 21 posterior mean (as Note 45b), driven by measured CO2 and light. Air T per growth sample:
#   meas       measured air T (reference)
#   n54        n54 closed-loop model air T, full season from day 0 (the coupled forecast situation)
#   n54 d9+    model T from day 9, measured T on days 0-8 (isolates the early-season failure of Note 56)
#   n54 OL     n54 with the logged vents replayed (contrast; Note 56 OL is warm by day)
# Outputs: final nodes, yield, yield d42-114 (Note 31 holdout), LAI, weekly-yield RMSE vs observations, and the
# node gap at days 8 / 48 / 100. This tells whether the crop parameters must be re-fitted (Option B) on n54.
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n54")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2")
    isfile(f) || error("missing $f")
    push!(parts, JLD2.load(f)["chain"][30001:end, 1:NM])
end
pm = vec(mean(reduce(vcat, parts), dims = 1)); pP = update_params(params_cl, CQ_MODEL, pm)
Tcl = full_eval(pm)[1]
usOL = simdays(update_params(pP, ["vent_ctrl"], [0.0]), days, u0_cl, RELTOL); usOL === nothing && error("OL solve failed")
Tol = usOL[2, :]

# Note 21 crop posterior mean, read directly from its chain (`centre` is redefined by the climate script)
c21 = JLD2.load(N21_CHAIN); m21 = vec(mean(c21["chain"][c21["burn_in"]:end, 1:end-1], dims = 1))
g21 = [m21[findfirst(==(q), c21["QoI"])] for q in QoI_G]
gp21 = update_params(gp_base, QoI_G, g21)
println("\nCrop parameters (Note 21 posterior mean): ", join(["$q=$(round(v, sigdigits=3))" for (q, v) in zip(QoI_G, g21)], "  "))

function run_T(src)          # src(d, k) -> air T [K] for growth sample k of day d (0-based)
    for d in 0:(N_DAYS - 1)
        M = work[d + 1]
        @inbounds for k in axes(M, 2); M[1, k] = src(d, k); end
    end
    simulate_growth_measured(params, gp21, season_B; LAI0 = 0.5)
end
meas(d, k)  = samp_meas[d + 1][1, k]
model(Tv)   = (d, k) -> (v = Tv[samp_col[d + 1][k]]; isnan(v) ? meas(d, k) : v)
yat(v, d)   = v[min(d, length(v) - 1) + 1]
function metrics(r)
    wk = weekly_incr(r.day, r.yield_FW, edges)
    (node = r.node[end], yld = r.yield_FW[end], yhold = yat(r.yield_FW, 114) - yat(r.yield_FW, 42), lai = r.LAI[end],
     wrmse = sqrt(mean(abs2, wk .- Δy_obs)), n8 = yat(r.node, 8), n48 = yat(r.node, 48), n100 = yat(r.node, 100))
end

# daily-mean T bias of each driver vs measured (growth samples)
function tbias(src)
    b = [mean(src(d, k) - meas(d, k) for k in axes(samp_meas[d + 1], 2)) for d in 0:(N_DAYS - 1) if size(samp_meas[d + 1], 2) > 0]
    mean(b)
end
runs = [("measured T", meas),
        ("n54 CL (coupled)", model(Tcl)),
        ("n54 CL, meas d0-8", (d, k) -> d < 9 ? meas(d, k) : model(Tcl)(d, k)),
        ("n54 OL (logged vents)", model(Tol))]
ref = nothing
println("\n", "="^150)
@printf("  %-24s | %7s | %7s %6s | %7s %6s | %7s %6s | %5s | %7s | %6s %6s %6s\n", "driver", "T bias", "nodes", "Δ%",
        "yield", "Δ%", "y42-114", "Δ%", "LAI", "wkRMSE", "Δn d8", "Δn d48", "Δn100")
println("  ", "-"^148)
pct(a, b) = 100 * (a - b) / b
for (nm_, src) in runs
    m = metrics(run_T(src))
    global ref; ref === nothing && (ref = m)
    @printf("  %-24s | %+7.2f | %7.1f %+6.1f | %7.2f %+6.1f | %7.2f %+6.1f | %5.2f | %7.3f | %+6.2f %+6.2f %+6.2f\n", nm_,
            tbias(src), m.node, pct(m.node, ref.node), m.yld, pct(m.yld, ref.yld), m.yhold, pct(m.yhold, ref.yhold),
            m.lai, m.wrmse, m.n8 - ref.n8, m.n48 - ref.n48, m.n100 - ref.n100)
end
@printf("\n  Observed final yield %.2f kg/m2 (Production.csv)\n", y_obs[end])
println("""
  T bias = mean daily air-T difference vs measured over growth samples (K). Δ% vs the measured-T run.
  Δn = node-number difference vs the measured-T run at days 8, 48 and 100.
  Read: if the coupled run stays within ~2-5 % nodes and the weekly-yield RMSE does not degrade, the Note 21 crop
  parameters can be kept with n54; otherwise re-fit them (Option B, long run -- ask first).""")
