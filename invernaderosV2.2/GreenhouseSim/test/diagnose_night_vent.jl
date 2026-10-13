# =============================================================================
# diagnose_night_vent.jl -- is the +0.47 K night warm bias a ventilation deficit or the soil anchor? (Note 38)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     julia --project=. test/diagnose_night_vent.jl
#
# Starts from the Note 38 pooled posterior mean (PRIORS=n37, closed-loop blocks, seeds 1+2). New default-off
# engine options (forward_map.jl, Note 38): vc_night_log (replay the logged vent at night) and vc_min_night
# (minimum night opening); night = Iglobal < vc_night_I (1 W/m2).
#   R0 posterior mean        closed-loop vents day and night (reproduces Note 38)
#   R1 hybrid                closed loop by day, LOGGED vents at night
#   R2 min night = median    closed loop + minimum night opening = median logged night opening (TRAIN days)
#   R3 min night = median/2
#   R4 soil lever            posterior mean with k_ground = 10 (vents unchanged)
# Reading: if R1/R2 remove the night bias without hurting the day, night warmth is a ventilation deficit
# (-> minimum-vent flag); if only R4 moves it, it is the soil anchor (-> soil ladder).
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n37"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
for k in (:vc_night_log, :vc_min_night, :vc_night_I)
    hasproperty(params, k) || error("constants_climate.json has no $k entry (Note 38 flags missing)")
end

# ── Note 38 pooled posterior mean ───────────────────────────────────────────
parts = Vector{Vector{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n37_seed$(s).jld2")
    isfile(f) || error("missing $f (Note 38 chain)")
    d = JLD2.load(f); push!(parts, vec(mean(d["chain"][d["burn_in"]+1:end, 1:NM], dims = 1)))
end
pm = (parts[1] .+ parts[2]) ./ 2
println("\nNote 38 pooled posterior mean: ", join(["$q=$(round(v, sigdigits=3))" for (q, v) in zip(CQ_MODEL, pm)], "  "))

# median logged night opening on TRAIN days only (no peeking at TEST)
rowday = [floor(Int, (t - t0) / 86400) for t in ts]
nightrows = [i for i in eachindex(ts) if Ig_v[i] < 1.0 && rowday[i] in TRAIN_DAYS && !isnan(U8log[i])]
U8med = median(U8log[nightrows]); U8mean = mean(U8log[nightrows])
@printf("Logged night vent opening on TRAIN days: median %.3f, mean %.3f  (%d rows)\n", U8med, U8mean, length(nightrows))

function run_p(names, vals)
    p = update_params(params_cl, CQ_MODEL, pm)
    isempty(names) || (p = update_params(p, names, vals))
    us = simdays(p, days, u0_8, RELTOL); us === nothing && error("solver failed for $(names)")
    us[2, :], rhf.(us[2, :], us[4, :]), clamp.(us[7, :], 0.0, 1.0)
end
ix(q) = findfirst(==(q), CQ_MODEL)

# check 0: flags explicitly 0 must be bit-for-bit identical to the JSON defaults
T0a, _, _ = run_p(String[], Float64[])
T0b, _, _ = run_p(["vc_night_log", "vc_min_night"], [0.0, 0.0])
@printf("[check 0] flags off (JSON) vs flags explicitly 0: max |ΔT| = %.3e K (%s)\n",
        maximum(abs.(T0a .- T0b)), T0a == T0b ? "bit-for-bit identical" : "NOT identical")

runs = [("R0 posterior",       String[], Float64[]),
        ("R1 night = logged",  ["vc_night_log"], [1.0]),
        ("R2 min night = med", ["vc_min_night"], [U8med]),
        ("R3 min night = med/2", ["vc_min_night"], [U8med / 2]),
        ("R4 k_ground = 10",   ["k_ground"], [10.0])]
res = Any[]
println("Running $(length(runs)) full-season simulations (reltol $(RELTOL)) ...")
for (nm_, nn, vv) in runs
    t = @elapsed (out = run_p(nn, vv)); push!(res, (nm_, out...)); @printf("  %-22s %.1f s\n", nm_, t)
end

# regression check against the Note 38 log (COMMON TEST 1.33 / +0.12 / +0.47 / -0.20)
eval_table(res[1][2], res[1][3], res[1][4]; label = "R0 regression check (should match the Note 38 log):")

# ── metrics ──────────────────────────────────────────────────────────────────
okT = .!isnan.(Tobs_g); okT[1] = false; okU = .!isnan.(U8obs); okU[1] = false
okV = okT .& .!isnan.(RHobs_g)
inset(ds) = in.(cd_, Ref(Set(ds)))
COM = inset(COMMON_TEST); EARLY = inset(filter(<(63), COMMON_TEST)); LATE = inset(filter(>=(63), COMMON_TEST))
mT = okT .& COM; mU = okU .& COM; mV = okV .& COM
VPo = RHobs_g ./ 100 .* Pws.(Tobs_g) ./ 100
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
rmse(v, o, m) = any(m) ? sqrt(mean((v[m] .- o[m]) .^ 2)) : NaN
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
function hourly_rmse(v, o, mask)
    e = Float64[]
    for cols in values(HB); k = [c for c in cols if mask[c]]; length(k) >= 6 && push!(e, mean(v[k]) - mean(o[k])); end
    isempty(e) ? NaN : sqrt(mean(abs2, e))
end
function wsd(v, mask)
    s = Float64[]
    for cols in values(HB); k = [c for c in cols if mask[c]]; length(k) >= 9 && push!(s, std(v[k])); end
    isempty(s) ? NaN : mean(s)
end

println("\n", "="^134)
println("  COMMON TEST. T in K, VP in hPa, vent 0-1. 'n early/late' = night T bias on common-test days < 63 / >= 63.")
@printf("  %-22s | %6s %6s %6s %6s %6s %7s %7s | %6s %6s | %7s %7s | %6s\n", "run", "T hRMS", "T 5m", "bias",
        "night", "day", "n early", "n late", "VPday", "VPnite", "U8 nite", "U8 day", "wsdN")
println("  ", "-"^132)
@printf("  %-22s | %6s %6s %6s %6s %6s %7s %7s | %6s %6s | %7s %7s | %6.3f\n", "observed", "", "", "", "", "", "", "",
        "", "", "", "", wsd(replace(Tobs_g, NaN => 0.0), mT .& nite))
for (nm_, T, RH, U) in res
    VPm = RH ./ 100 .* Pws.(T) ./ 100
    @printf("  %-22s | %6.3f %6.3f %+6.2f %+6.2f %+6.2f %+7.2f %+7.2f | %+6.2f %+6.2f | %+7.3f %+7.3f | %6.3f\n", nm_,
            hourly_rmse(T, Tobs_g, mT), rmse(T, Tobs_g, mT), mb(T, Tobs_g, mT), mb(T, Tobs_g, mT .& nite),
            mb(T, Tobs_g, mT .& dayt), mb(T, Tobs_g, okT .& EARLY .& nite), mb(T, Tobs_g, okT .& LATE .& nite),
            mb(VPm, VPo, mV .& dayt), mb(VPm, VPo, mV .& nite),
            mb(U, U8obs, mU .& nite), mb(U, U8obs, mU .& dayt), wsd(T, mT .& nite))
end
println("""
  How to read it:
   - R1/R2 night bias -> ~0 with day unchanged  : night warmth = ventilation deficit (minimum-vent flag next).
   - only R4 moves the night bias               : soil anchor (soil-ladder flag next).
   - both move it                               : the two are confounded; calibrate with the vent flag first
                                                  (it is directly observed through U8 night bias).
   VP night tells whether the extra night ventilation also fixes the humid nights (+0.43 hPa).""")

# ── figure: four nights ──────────────────────────────────────────────────────
SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "56,73,93,105"), ","))
ps = Any[]
for d in SHOWD
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    pT = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 3, label = "measured", title = "day $d", ylabel = "T air (°C)",
              legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6)
    pV = plot(h, U8obs[c], lc = :black, lw = 2, label = "logged", xlabel = "hour", ylabel = "U8", ylims = (-0.05, 0.8), legend = false)
    for (k, (nm_, T, _, U)) in enumerate(res)
        k in (1, 2, 3, 5) || continue
        plot!(pT, h, T[c] .- 273.15, lw = 1.2, label = nm_); plot!(pV, h, U[c], lw = 1.2, label = nm_)
    end
    push!(ps, pT); push!(ps, pV)
end
order = [ps[2(i-1)+r] for r in 1:2 for i in 1:length(SHOWD)]
fig = joinpath(FIGDIR, "AiCU_night_vent_diagnostic.png")
savefig(plot(order..., layout = (2, length(SHOWD)), size = (420 * length(SHOWD), 650)), fig)
println("Figure: ", fig)
