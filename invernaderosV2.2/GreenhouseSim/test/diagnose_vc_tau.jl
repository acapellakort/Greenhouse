# =============================================================================
# diagnose_vc_tau.jl -- is vc_tau -> 2 h a likelihood blind spot or a thermal-inertia compensation?
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1-2 min, no MCMC):
#     julia --project=. test/diagnose_vc_tau.jl
# Part A: tau sweep.  Part B: post-hoc sensor lag on model T.  Part C: limit-cycle check (gain, vc_hband).
#
# Starts from the Note 36 pooled posterior mean (closed-loop blocks, seeds 1+2) and changes ONLY the
# vent-controller time constant (and, in V5, the controller to the Note 34 identified values):
#   V1 posterior mean  (vc_tau ~ 6960 s)
#   V2 vc_tau 1800 s
#   V3 vc_tau  900 s
#   V4 vc_tau  300 s
#   V5 identified ctrl: vc_tau 300 s, vc_gain 0.11/K, vc_e0 -2.5 K (Note 34)
# Full season, closed loop, reltol 1e-4.
#
# Reading the result:
#   * 5-min vent RMSE and within-hour vent variability should favour SHORT tau if the grower's vents
#     really react in <= 5 min (Note 34).
#   * If T (hourly RMSE, within-hour T variability, spike overshoot) favours LONG tau at the same time,
#     the slow vent is compensating missing thermal inertia (too-fast air temperature), not identifying
#     the controller.
#   * If T is indifferent to tau, the hourly likelihood is simply blind to tau.
# =============================================================================

ENV["CL_MODE"] = "setup"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

# ── Note 36 pooled posterior mean ───────────────────────────────────────────
pm_parts = Vector{Vector{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_seed$(s).jld2")
    isfile(f) || error("missing $f (Note 36 chain)")
    d = JLD2.load(f); ch = d["chain"]; b = d["burn_in"]
    push!(pm_parts, vec(mean(ch[b+1:end, 1:NM], dims = 1)))
end
pm = (pm_parts[1] .+ pm_parts[2]) ./ 2
println("\nNote 36 pooled posterior mean: ", join(["$q=$(round(v, sigdigits=3))" for (q, v) in zip(CQ_MODEL, pm)], "  "))

ix(q) = findfirst(==(q), CQ_MODEL)
withv(x, pairs...) = (y = copy(x); for (q, v) in pairs; y[ix(q)] = v; end; y)
variants = [("V1 post. mean",     pm),
            ("V2 tau 1800",       withv(pm, "vc_tau" => 1800.0)),
            ("V3 tau 900",        withv(pm, "vc_tau" => 900.0)),
            ("V4 tau 300",        withv(pm, "vc_tau" => 300.0)),
            ("V5 identified ctl", withv(pm, "vc_tau" => 300.0, "vc_gain" => 0.11, "vc_e0" => -2.5))]

res = Any[]
println("Running $(length(variants)) full-season simulations (reltol $(RELTOL)) ...")
for (nm_, x) in variants
    t = @elapsed (us = simdays(update_params(params_cl, CQ_MODEL, x), days, u0_8, RELTOL))
    us === nothing && error("$nm_: solver failed")
    push!(res, (nm_, us[2, :], clamp.(us[7, :], 0.0, 1.0)))
    @printf("  %-18s %.1f s\n", nm_, t)
end

# ── metrics ──────────────────────────────────────────────────────────────────
okT = .!isnan.(Tobs_g); okT[1] = false
okU = .!isnan.(U8obs);  okU[1] = false
inset(ds) = in.(cd_, Ref(Set(ds)))
COM = inset(COMMON_TEST)
hr(c) = (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600))
# group columns into (day, hour) bins once
const HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, hr(c), Int[]), c); end

# mean within-hour standard deviation of v over bins whose columns are all in `mask` (>= 9 of 12 valid)
function within_hour_sd(v, mask)
    s = Float64[]
    for cols in values(HB)
        k = [c for c in cols if mask[c]]; length(k) >= 9 || continue
        push!(s, std(v[k]))
    end
    isempty(s) ? NaN : mean(s)
end
# 1-hour increment SD (12 steps of 5 min, same day)
function incr_sd(v, mask)
    d = Float64[]
    for c in 2:NCOL-12
        (mask[c] && mask[c+12] && cd_[c] == cd_[c+12]) || continue
        push!(d, v[c+12] - v[c])
    end
    isempty(d) ? NaN : std(d)
end
function hourly_rmse(vm, vo, mask)
    e = Float64[]
    for cols in values(HB)
        k = [c for c in cols if mask[c]]; length(k) >= 6 || continue
        push!(e, mean(vm[k]) - mean(vo[k]))
    end
    isempty(e) ? NaN : sqrt(mean(abs2, e))
end
rmse(v, o, m) = any(m) ? sqrt(mean((v[m] .- o[m]) .^ 2)) : NaN
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
function overshoot(T)
    ov = Float64[]
    for d in 11:41
        m = okT .& dayt .& (cd_ .== d); count(m) < 12 && continue
        push!(ov, maximum(T[m]) - maximum(Tobs_g[m]))
    end
    mean(ov), count(>(3), ov), length(ov)
end

mT = okT .& COM; mTd = mT .& dayt; mU = okU .& COM; mUd = mU .& dayt
Tobs_z = replace(Tobs_g, NaN => 0.0); Uobs_z = replace(U8obs, NaN => 0.0)
obs_wsdT = within_hour_sd(Tobs_z, mTd); obs_incT = incr_sd(Tobs_z, mTd)
obs_wsdU = within_hour_sd(Uobs_z, mUd); obs_incU = incr_sd(Uobs_z, mUd)

println("\n", "="^128)
println("  COMMON TEST days (49-62, 70-76, 84-90, 98-114). T in K, vent opening 0-1. 'wsd' = mean within-hour SD (daytime),")
println("  'Δ1h' = SD of 1-hour increments (daytime). Observed values in the first row. Spike overshoot: days 11-41, daytime max.")
@printf("  %-18s | %6s %6s %6s %6s %6s %6s | %6s %6s %6s %6s %6s | %6s %5s\n", "", "T hRMS", "T 5m", "night", "day",
        "T wsd", "T Δ1h", "U 5m", "U hRMS", "U bias", "U wsd", "U Δ1h", "ovrsh", ">3K")
println("  ", "-"^126)
@printf("  %-18s | %6s %6s %6s %6s %6.3f %6.3f | %6s %6s %6s %6.3f %6.3f | %6s %5s\n", "observed", "", "", "", "",
        obs_wsdT, obs_incT, "", "", "", obs_wsdU, obs_incU, "", "")
for (nm_, T, U) in res
    ov, n3, nd = overshoot(T)
    @printf("  %-18s | %6.3f %6.3f %+6.2f %+6.2f %6.3f %6.3f | %6.3f %6.3f %+6.3f %6.3f %6.3f | %+6.2f %2d/%2d\n", nm_,
            hourly_rmse(T, Tobs_z, mT), rmse(T, Tobs_g, mT), mb(T, Tobs_g, mT .& nite), mb(T, Tobs_g, mTd),
            within_hour_sd(T, mTd), incr_sd(T, mTd),
            rmse(U, U8obs, mU), hourly_rmse(U, Uobs_z, mU), mb(U, U8obs, mU),
            within_hour_sd(U, mUd), incr_sd(U, mUd), ov, n3, nd)
end
println("""
  How to read it:
   - U 5m / U wsd closer to observed for short tau  -> the 5-min data support a fast controller (Note 34).
   - T columns better (or T wsd closer to observed) for long tau -> vc_tau compensates missing thermal inertia.
   - T columns ~unchanged across tau                -> the hourly likelihood is blind to vc_tau.
   V2-V5 keep all other parameters at the V1 posterior mean (tuned for tau ~2 h), so they are a sensitivity
   test, not a re-calibration.""")

# ── figure: T and vent on four days ──────────────────────────────────────────
SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "24,30,63,93"), ","))
ps = Any[]
for d in SHOWD
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    pT = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 3, label = "measured", title = "day $d", ylabel = "T air (°C)",
              legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6)
    pV = plot(h, U8obs[c], lc = :black, lw = 2, label = "logged vent", xlabel = "hour", ylabel = "U8", ylims = (-0.05, 1.05),
              legend = false)
    for (k, (nm_, T, U)) in enumerate(res)
        k in (1, 3, 4, 5) || continue
        plot!(pT, h, T[c] .- 273.15, lw = 1.3, label = nm_); plot!(pV, h, U[c], lw = 1.3, label = nm_)
    end
    push!(ps, pT); push!(ps, pV)
end
order = [ps[2(i-1)+r] for r in 1:2 for i in 1:length(SHOWD)]
fig = joinpath(FIGDIR, "AiCU_vc_tau_diagnostic.png")
savefig(plot(order..., layout = (2, length(SHOWD)), size = (420 * length(SHOWD), 650)), fig)
println("Figure: ", fig)

# =============================================================================
# PART B -- observation operator: first-order sensor lag applied to the model T (post hoc)
# =============================================================================
# y_k = y_{k-1} + a (x_k - y_{k-1}),  a = 1 - exp(-300 s / tau_s), on the 5-min grid; the duplicated
# day-boundary column (j = 0 of each day after the first) is a zero-length step (a = 0).
function sensor_lag(x, tau_s)
    tau_s <= 0 && return copy(x)
    a = 1 - exp(-300.0 / tau_s); y = similar(x); y[1] = x[1]; y[2] = x[2]
    for c in 3:length(x)
        j = (c - 2) % NPD
        y[c] = j == 0 ? y[c-1] : y[c-1] + a * (x[c] - y[c-1])
    end
    y
end
mTn = mT .& nite; mUn = mU .& nite
obs_wsdTn = within_hour_sd(Tobs_z, mTn)
println("\n", "="^110)
println("  PART B -- model T passed through a first-order SENSOR LAG tau_s (post hoc). COMMON TEST.")
println("  If tau_s ~ 5-10 min brings 'T wsd' and 'T Δ1h' to the observed values without hurting hRMS, the excess")
println("  within-hour variability is a measurement effect (smoothing sensor / averaging), not missing physics.")
@printf("  %-18s %6s | %6s %6s %6s %6s %7s\n", "run", "tau_s", "T hRMS", "T 5m", "T wsd", "T Δ1h", "wsd nite")
println("  ", "-"^70)
@printf("  %-18s %6s | %6s %6s %6.3f %6.3f %7.3f\n", "observed", "", "", "", obs_wsdT, obs_incT, obs_wsdTn)
for k in (1, 2, 3, 5)
    nm_, T, _ = res[k]
    for ts_ in (0.0, 300.0, 600.0, 1200.0)
        Tf = sensor_lag(T, ts_)
        @printf("  %-18s %4.0f m | %6.3f %6.3f %6.3f %6.3f %7.3f\n", ts_ == 0 ? nm_ : "", ts_ / 60,
                hourly_rmse(Tf, Tobs_z, mT), rmse(Tf, Tobs_g, mT), within_hour_sd(Tf, mTd), incr_sd(Tf, mTd),
                within_hour_sd(Tf, mTn))
    end
end

# =============================================================================
# PART C -- limit-cycle check: softer humidity band and lower vent gain
# =============================================================================
#   C1 tau 900, gain 0.15                       (gain only)
#   C2 tau 900, gain 0.15, vc_hband 3 g/m3      (gain + soft humidity band)
#   C3 tau 900, posterior gain, vc_hband 3      (band only)
#   C4 identified ctl (tau 300, 0.11, -2.5), vc_hband 3
# Night within-hour SD of T and vent is the limit-cycle indicator (observed values in the first row).
function simx(x; names = String[], vals = Float64[])
    p = update_params(params_cl, CQ_MODEL, x)
    isempty(names) || (p = update_params(p, names, vals))
    simdays(p, days, u0_8, RELTOL)
end
x900 = withv(pm, "vc_tau" => 900.0)
cvars = [("V3 tau 900 (ref)",  x900, String[], Float64[]),
         ("C1 g0.15",          withv(x900, "vc_gain" => 0.15), String[], Float64[]),
         ("C2 g0.15 hb3",      withv(x900, "vc_gain" => 0.15), ["vc_hband"], [3.0]),
         ("C3 hb3",            x900, ["vc_hband"], [3.0]),
         ("C4 ident. hb3",     withv(pm, "vc_tau" => 300.0, "vc_gain" => 0.11, "vc_e0" => -2.5), ["vc_hband"], [3.0])]
resC = Any[]
println("\nRunning $(length(cvars)) limit-cycle variants ...")
for (nm_, x, nn, vv) in cvars
    us = simx(x; names = nn, vals = vv); us === nothing && error("$nm_: solver failed")
    push!(resC, (nm_, us[2, :], clamp.(us[7, :], 0.0, 1.0)))
end
obs_wsdUn = within_hour_sd(Uobs_z, mUn)
println("\n", "="^120)
println("  PART C -- limit-cycle check. COMMON TEST. 'nite' columns = within-hour SD at night (the sawtooth indicator).")
@printf("  %-18s | %6s %6s %6s %6s %6s %7s | %6s %6s %6s %7s | %6s %5s\n", "run", "T hRMS", "T 5m", "night", "day",
        "T wsd", "T wsdN", "U 5m", "U bias", "U wsd", "U wsdN", "ovrsh", ">3K")
println("  ", "-"^118)
@printf("  %-18s | %6s %6s %6s %6s %6.3f %7.3f | %6s %6s %6.3f %7.3f | %6s %5s\n", "observed", "", "", "", "",
        obs_wsdT, obs_wsdTn, "", "", obs_wsdU, obs_wsdUn, "", "")
for (nm_, T, U) in resC
    ov, n3, nd = overshoot(T)
    @printf("  %-18s | %6.3f %6.3f %+6.2f %+6.2f %6.3f %7.3f | %6.3f %+6.3f %6.3f %7.3f | %+6.2f %2d/%2d\n", nm_,
            hourly_rmse(T, Tobs_z, mT), rmse(T, Tobs_g, mT), mb(T, Tobs_g, mTn), mb(T, Tobs_g, mTd),
            within_hour_sd(T, mTd), within_hour_sd(T, mTn), rmse(U, U8obs, mU), mb(U, U8obs, mU),
            within_hour_sd(U, mUd), within_hour_sd(U, mUn), ov, n3, nd)
end
println("""
  How to read it: if C2/C3/C4 bring 'T wsdN' and 'U wsdN' down to the observed values (night sawtooth gone),
  the cycle comes from the controller shape (sharp humidity band / high gain), not from missing thermal inertia.""")

SHOWC = parse.(Int, split(get(ENV, "SHOW_DAYS_C", "93,100"), ","))
psC = Any[]
for d in SHOWC
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    pT = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 3, label = "measured", title = "day $d (Part C)",
              ylabel = "T air (°C)", legend = d == SHOWC[1] ? :topleft : false, legendfontsize = 6)
    pV = plot(h, U8obs[c], lc = :black, lw = 2, label = "logged", xlabel = "hour", ylabel = "U8", ylims = (-0.05, 1.05), legend = false)
    for (nm_, T, U) in resC; plot!(pT, h, T[c] .- 273.15, lw = 1.2, label = nm_); plot!(pV, h, U[c], lw = 1.2, label = nm_); end
    push!(psC, pT); push!(psC, pV)
end
orderC = [psC[2(i-1)+r] for r in 1:2 for i in 1:length(SHOWC)]
figC = joinpath(FIGDIR, "AiCU_vc_tau_limitcycle.png")
savefig(plot(orderC..., layout = (2, length(SHOWC)), size = (520 * length(SHOWC), 650)), figC)
println("Figure: ", figC)
