# =============================================================================
# postprocess_joint_optionB.jl -- posteriors + full-season (0..114 d) prediction
# =============================================================================
# Run from the GreenhouseSim folder (after the `full` runs have finished):
#     julia --project=. test/postprocess_joint_optionB.jl
#     NDRAW=100 julia --project=. test/postprocess_joint_optionB.jl     # more posterior draws
#
# 1. Loads every AiCU_joint_optionB_d42_seed*.jld2 chain (post burn-in), plots the
#    marginal posteriors with seeds overlaid, prints split-R-hat and flags parameters
#    whose posterior sits on a prior bound.
# 2. Re-uses the joint script's input builder for the FULL season (OPTB_MODE=setup,
#    N_DAYS=115) and runs Vanthoor -> growth for NDRAW posterior draws + the mean.
#    Days 0-41 = calibration window; days 42-114 = HOLDOUT (never seen by the chain;
#    outdoor weather and grower setpoints are the observed ones).
# 3. Compares three ways to predict the season:
#      B  joint posterior (Vanthoor T drives growth)              -- the Option B result
#      S  Note 21 growth params + Vanthoor T at joint theta_v mean  -- "sequential"
#      A  Note 21 growth params + MEASURED T                        -- in-sample for growth!
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = get(ENV, "N_FULL", "115")
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))   # full-season inputs, forward map, growth obs

const CAL_DAYS = 42
const NDRAW    = parse(Int, get(ENV, "NDRAW", "40"))
const PP_TAG   = "AiCU_joint_optionB_d42_seed"
const PP_FIG   = joinpath(FIGDIR, "AiCU_joint_optionB_pp")

# ── 1. chains ────────────────────────────────────────────────────────────────
pp_files = sort(filter(f -> startswith(f, PP_TAG) && endswith(f, ".jld2") && !occursin("ckpt", f), readdir(POSTDIR)))
isempty(pp_files) && error("no joint chains found in $POSTDIR")
pp_ch = [JLD2.load(joinpath(POSTDIR, f)) for f in pp_files]
all(c["QoI"] == QoI for c in pp_ch) || error("chain QoI does not match the script's QoI")
posts = [c["chain"][c["burn_in"]+1:end, 1:end-1] for c in pp_ch]
println("\nChains: ", join(pp_files, ", "), "   (", size(posts[1], 1), " post-burn-in rows each)")

# split R-hat: each chain cut in halves -> 2m sequences
function split_rhat(cols)
    seqs = Vector{Vector{Float64}}()
    for c in cols
        h = length(c) ÷ 2; push!(seqs, c[1:h]); push!(seqs, c[h+1:2h])
    end
    n = length(seqs[1]); W = mean(var.(seqs)); B = n * var(mean.(seqs))
    return sqrt(((n - 1) / n * W + B / n) / W)
end

pool = vcat(posts...)
θmean = vec(mean(pool, dims = 1))
@printf("\n%-13s %11s %11s %11s   %6s   %s\n", "param", "mean", "2.5%", "97.5%", "Rhat", "flag")
hists = Any[]
for (i, q) in enumerate(QoI)
    lo, _, hi = QoI_dict[q]
    qs = quantile(pool[:, i], [0.025, 0.975])
    rh = split_rhat([p[:, i] for p in posts])
    tr = (lo > 0 && hi / lo > 100) ? log : identity          # wide positive priors (nu4) judged in log space
    span = tr(hi) - tr(lo)
    flag = tr(qs[1]) - tr(lo) < 0.02span ? "AT LOWER BOUND" : tr(hi) - tr(qs[2]) < 0.02span ? "AT UPPER BOUND" : ""
    rh > 1.1 && (flag *= (isempty(flag) ? "" : "; ") * "not mixed")
    @printf("%-13s %11.4g %11.4g %11.4g   %6.3f   %s\n", q, θmean[i], qs[1], qs[2], rh, flag)
    h = plot(title = @sprintf("%s  (R̂=%.2f)", q, rh), titlefontsize = 9, legend = i == 1 ? :topright : false)
    for (k, p) in enumerate(posts)
        histogram!(h, p[:, i], bins = 40, normalize = :pdf, alpha = 0.45, label = "seed $(pp_ch[k]["seed"])", lw = 0)
    end
    vline!(h, [centre[i]], lc = :black, ls = :dash, label = "warm start")
    !isempty(flag) && startswith(flag, "AT") && vline!(h, [startswith(flag, "AT LOWER") ? lo : hi], lc = :red, lw = 2, label = "prior bound")
    push!(hists, h)
end
savefig(plot(hists..., layout = length(hists), size = (1500, 1150)), PP_FIG * "_posteriors.png")

# ── 2. full-season prediction ────────────────────────────────────────────────
function run_full(θv, θg)
    T2 = try forward_map(θv, QoI_V, base, days, tspan, teval; u0 = copy(u0))[2] catch; return nothing end
    T2 = collect(T2)
    (length(T2) == NCOL && all(isfinite, T2)) || return nothing
    for d in 1:N_DAYS
        M = work[d]; c = samp_col[d]
        @inbounds for k in eachindex(c); M[1, k] = T2[c[k]]; end
    end
    res = simulate_growth_measured(params, update_params(gp_base, QoI_G, θg), season_B; LAI0 = 0.5)
    return (T2 = T2, yield = copy(res.yield_FW), node = copy(res.node), LAI = copy(res.LAI))
end

rng = MersenneTwister(2026)
draws = pool[rand(rng, 1:size(pool, 1), NDRAW), :]
println("\nSimulating full season ($(N_DAYS) days) for posterior mean + $NDRAW draws ...")
t_pp = @elapsed begin
    simB  = run_full(θmean[IV], θmean[IG])
    simBd = filter(!isnothing, [run_full(draws[k, IV], draws[k, IG]) for k in 1:NDRAW])
    simS  = run_full(θmean[IV], centre[IG])
end
@printf("  done in %.0f s (%d/%d draws finite)\n", t_pp, length(simBd), NDRAW)
resA = simulate_growth_measured(params, update_params(gp_base, QoI_G, centre[IG]), season_A; LAI0 = 0.5)
simA = (T2 = nothing, yield = resA.yield_FW, node = resA.node, LAI = resA.LAI)

# ── 3. metrics: calibration (0-41) vs holdout (42-end) ───────────────────────
coldays = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
okT  = .!isnan.(Tobs_g); okT[1] = false
mcal = okT .& (coldays .< CAL_DAYS); mhol = okT .& (coldays .>= CAL_DAYS)
rmse(T, m) = sqrt(mean((T[m] .- Tobs_g[m]) .^ 2)); bias(T, m) = mean(T[m] .- Tobs_g[m])
obs_cum(d) = cum_at(days_obs, y_obs, d)
hol_edges = filter(>=(CAL_DAYS - 1), edges)
wk_rmse(y) = sqrt(mean((weekly_incr(0:N_DAYS-1, y, hol_edges) .- weekly_incr(days_obs, y_obs, hol_edges)) .^ 2))
hol_cm = cm_days .>= CAL_DAYS
node_rmse(nd) = sqrt(mean((cm_N[hol_cm] .- [at_day(nd, d) for d in cm_days[hol_cm]]) .^ 2))

println("\n", "="^100)
@printf("  AiCU full season: calibration days 0-%d, HOLDOUT days %d-%d\n", CAL_DAYS - 1, CAL_DAYS, LAST_DAY)
println("="^100)
@printf("  %-34s %8s %8s %8s | %9s %10s %9s | %9s\n", "", "RMSE_T", "RMSE_T", "bias_T", "yield", "holdout", "wk-yield", "nodes")
@printf("  %-34s %8s %8s %8s | %9s %10s %9s | %9s\n", "", "cal [K]", "hold [K]", "hold [K]", "d$(LAST_DAY)", "Δyield", "RMSE hold", "RMSE hold")
@printf("  %-34s %8s %8s %8s | %9.2f %10.2f %9s | %9s\n", "observed", "", "", "", obs_cum(LAST_DAY), obs_cum(LAST_DAY) - obs_cum(CAL_DAYS - 1), "", "")
function row(name, s)
    Tc = s.T2 === nothing ? "--" : @sprintf("%.2f", rmse(s.T2, mcal))
    Th = s.T2 === nothing ? "--" : @sprintf("%.2f", rmse(s.T2, mhol))
    Tb = s.T2 === nothing ? "--" : @sprintf("%+.2f", bias(s.T2, mhol))
    @printf("  %-34s %8s %8s %8s | %9.2f %10.2f %9.3f | %9.2f\n", name, Tc, Th, Tb,
            s.yield[LAST_DAY + 1], s.yield[LAST_DAY + 1] - s.yield[CAL_DAYS], wk_rmse(s.yield), node_rmse(s.node))
end
row("B  joint posterior mean", simB)
row("S  N21 growth + Vanthoor T", simS)
row("A  N21 growth + measured T [in-sample]", simA)
Yd = reduce(hcat, [s.yield for s in simBd])
q5, q95 = quantile(Yd[LAST_DAY + 1, :], [0.05, 0.95])
@printf("\n  B draws: final yield 90%% band [%.2f, %.2f] kg/m2 (parameter uncertainty only, no obs noise)\n", q5, q95)
@printf("  Observed final yield %s the B band.\n", q5 <= obs_cum(LAST_DAY) <= q95 ? "is INSIDE" : "is OUTSIDE")

# ── 4. figures (x axis in days) ──────────────────────────────────────────────
dayx = 0:N_DAYS-1
_q(v, p) = any(isnan, v) ? NaN : quantile(v, p)
band(M) = ([_q(M[i, :], 0.05) for i in axes(M, 1)], [_q(M[i, :], 0.95) for i in axes(M, 1)])

# daily-mean temperature (5-min traces are unreadable over 115 days)
dmean(v) = [let m = okT .& (coldays .== d); any(m) ? mean(v[m]) - 273.15 : NaN end for d in dayx]
TdB = reduce(hcat, [dmean(s.T2) for s in simBd]); lo_, hi_ = band(TdB)
pT = plot(dayx, dmean(Tobs_g), lc = :black, lw = 2, label = "measured", xlabel = "day", ylabel = "daily mean T air (°C)",
          title = "Daily mean T_air — full season", size = (1000, 380), legend = :bottomleft)
plot!(pT, dayx, lo_, fillrange = hi_, fillalpha = 0.25, lw = 0, c = :darkorange, label = "B 90% band")
plot!(pT, dayx, dmean(simB.T2), lc = :darkorange, lw = 2, label = "B posterior mean")
vspan!(pT, [CAL_DAYS - 0.5, N_DAYS], alpha = 0.06, c = :blue, label = "holdout")
savefig(pT, PP_FIG * "_T_daily.png")

dx5 = (collect(1:NCOL) .- 2) ./ NPD
pT5 = plot(dx5, Tobs_g .- 273.15, lc = :black, lw = 1, label = "measured", xlabel = "day", ylabel = "T air (°C)",
           title = "T_air 5-min — full season", size = (1600, 400))
plot!(pT5, dx5, simB.T2 .- 273.15, lc = :darkorange, lw = 0.7, label = "B posterior mean")
vline!(pT5, [CAL_DAYS], lc = :blue, ls = :dash, label = "end calibration")
savefig(pT5, PP_FIG * "_T_5min.png")

lo_, hi_ = band(Yd)
pY = plot(days_obs, y_obs, lc = :black, lw = 3, label = "observed", xlabel = "day", ylabel = "cumulative yield (kg FW/m²)",
          title = "Cumulative yield — full season", legend = :topleft, size = (900, 450))
plot!(pY, dayx, lo_, fillrange = hi_, fillalpha = 0.25, lw = 0, c = :darkorange, label = "B 90% band")
plot!(pY, dayx, simB.yield, lc = :darkorange, lw = 2, label = "B joint")
plot!(pY, dayx, simS.yield, lc = :purple, lw = 2, ls = :dot, label = "S N21 + Vanthoor T")
plot!(pY, dayx, simA.yield, lc = :green, lw = 2, ls = :dash, label = "A N21 + measured T (in-sample)")
vline!(pY, [CAL_DAYS - 0.5], lc = :blue, ls = :dash, label = "end calibration")
savefig(pY, PP_FIG * "_yield.png")

Nd = reduce(hcat, [s.node for s in simBd]); lo_, hi_ = band(Nd)
pN = scatter(cm_days, cm_N, mc = :black, label = "observed", xlabel = "day", ylabel = "leaves / stem", title = "Node number — full season",
             legend = :topleft, size = (900, 450))
plot!(pN, dayx, lo_, fillrange = hi_, fillalpha = 0.25, lw = 0, c = :darkorange, label = "B 90% band")
plot!(pN, dayx, simB.node, lc = :darkorange, lw = 2, label = "B joint")
plot!(pN, dayx, simS.node, lc = :purple, lw = 2, ls = :dot, label = "S")
plot!(pN, dayx, simA.node, lc = :green, lw = 2, ls = :dash, label = "A (in-sample)")
vline!(pN, [CAL_DAYS - 0.5], lc = :blue, ls = :dash, label = "")
savefig(pN, PP_FIG * "_nodes.png")

println("\nFigures: ", PP_FIG, "_{posteriors,T_daily,T_5min,yield,nodes}.png")
