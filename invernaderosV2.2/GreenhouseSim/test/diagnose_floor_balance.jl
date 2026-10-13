# =============================================================================
# diagnose_floor_balance.jl -- floor energy bookkeeping, solar partition and floor-deep conductance (Note 46)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~3-5 min, no MCMC):
#     julia --project=. test/diagnose_floor_balance.jl
#
# Base structure = Note 44 run D3: Note 38 posterior mean, closed-loop vents, hc_log_pipe = 1, vc_min_log = 1
# (logged VentLeewMin input), pipe length gamma1 x 2, soil_dyn with sd_tau = 14 d.
# New default-off flags (Note 46):
#   floor_full  : r11 (floor FIR -> sky), r14 (canopy FIR -> floor, on T1) and h11 (air-floor conduction) are
#                 booked into the floor state instead of only into air / canopy
#   solar_split : the cover-absorbed sun (cov_abs * Iglob) heats the cover node, not the air (V2.0 puts
#                 1.0 * Iglob into the air via alpha8 = alpha9 = 1); the floor absorbs the transmitted sun
#                 that the canopy does not intercept
#   deep_k      : floor -> deep soil uses k_deep (2 W m-2 K-1) instead of k_ground
# Runs
#   R0  D3 reference (should reproduce Note 44 D3: hourly RMSE 1.254 K)
#   F1  + floor_full          K1  + deep_k          FK  + floor_full + deep_k
#   S1  + solar_split         SF  + solar_split + floor_full + deep_k
#   SFv SF with the vent gain halved (does the vent gain come back to the logged value?)
# Then a k_ground scan (3, 5, 10, 20, posterior) for R0, FK and SF: where is the best k_ground now?
# Parameters were calibrated with the OLD structure, so S-runs are expected to be off by day; they are a
# sensitivity test that decides which structure goes into the next calibration, not a result.
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n37"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
for k in (:hc_log_pipe, :vc_min_log, :soil_dyn, :sd_tau, :floor_full, :solar_split, :cov_abs, :deep_k, :k_deep)
    hasproperty(params, k) || error("constants_climate.json has no $k entry (Note 42-46 flags missing)")
end
isdefined(GreenhouseSim, :floor_full_fluxes) || error("src/climate.jl has no floor_full_fluxes (Note 46 edit missing)")

parts = Vector{Vector{Float64}}()
for s in (1, 2)
    d = JLD2.load(joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n37_seed$(s).jld2"))
    push!(parts, vec(mean(d["chain"][d["burn_in"]+1:end, 1:NM], dims = 1)))
end
pm = (parts[1] .+ parts[2]) ./ 2
p38 = update_params(params_cl, CQ_MODEL, pm)

VentMinv = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    vm = filter(!isnan, [_num(vip[i, Symbol("VentLeewMin_Vip_$z")]) for z in 1:6])
    push!(VentMinv, isempty(vm) ? 0.0 : clamp(mean(vm) / 200, 0.0, 1.0))
end
length(VentMinv) == length(ts) || error("row filter mismatch (VentMin)")
ctrl_v = ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                        Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv))
day0 = [i for i in eachindex(ts) if ts[i] < t0 + 86400 && !isnan(Tair_obs[i])]
u0_9 = vcat(u0_8, mean(Tair_obs[day0]))

function simrun(p, w, ctrl, u0x)
    us = Matrix{Float64}(undef, length(u0x), 1 + length(days) * NPD); us[:, 1] = u0x
    u = copy(u0x); colx = 2
    for d in days
        ctx = GreenhouseSim.SimContext(p, w, ctrl, t0 + 86400.0 * d, 86400.0 * d)
        sol = GreenhouseSim.solve(GreenhouseSim.ODEProblem(GreenhouseSim.rhs!, u, tspan, ctx), GreenhouseSim.Tsit5();
                                  saveat = teval, reltol = RELTOL)
        length(sol.t) == NPD || error("solver stopped on day $d")
        S = Array(sol); us[:, colx:colx+NPD-1] = S; colx += NPD; u = S[:, end]
    end
    us
end
pD3 = update_params(p38, ["hc_log_pipe", "vc_min_log", "gamma1", "soil_dyn", "sd_tau"],
                    [1.0, 1.0, 2.0 * p38.gamma1, 1.0, 14 * 86400.0])
fl(p; ff = 0.0, ss = 0.0, dk = 0.0, kg = p.k_ground) =
    update_params(p, ["floor_full", "solar_split", "deep_k", "k_ground"], [ff, ss, dk, kg])
@printf("\nNote 38 posterior: k_ground = %.1f, vc_gain = %.2f; fixed: k_deep = %.1f, cov_abs = %.2f, C_soil = %.2g\n",
        p38.k_ground, p38.vc_gain, p38.k_deep, p38.cov_abs, p38.C_soil)

# ── solar partition from the inputs alone (daylight means, W m^-2) ───────────
let P = p38, dl = [i for i in eachindex(Ig_v) if Ig_v[i] > 5.0]
    I2 = Ig_v[dl]; L = LAI_row[dl]
    old_air = I2 .* (P.eta1 * (P.tau1 * P.eta2 + (P.alpha2 + P.alpha7) * P.eta3) + (P.alpha8 * P.eta2 + P.alpha9 * P.eta3))
    new_air = I2 .* (P.eta1 * (P.tau1 * P.eta2 + (P.alpha2 + P.alpha7) * P.eta3))
    can_par = (1 - P.eta1) * P.tau1 * P.eta2 .* I2 .* ((1 - P.rho1) .* (1 .- exp.(-P.beta1 .* L)) .+
              exp.(-P.beta1 .* L) .* P.rho2 .* (1 - P.rho1) .* (1 .- exp.(-P.beta2 .* L)))
    can_nir = (1 - P.eta1) * P.alpha2 * P.eta3 .* I2
    flr     = (1 - P.eta1) * P.tau1 .* I2 .* (P.eta2 .* exp.(-P.beta1 .* L) .* (1 - P.rho2) .+ P.eta3 .* exp.(-P.beta3 .* L) .* P.alpha7)
    cov     = P.cov_abs .* I2
    println("\nSolar partition, daylight mean (W m^-2), outside global radiation mean = ", round(mean(I2), digits = 1))
    @printf("  %-34s %8s %8s\n", "", "V2.2", "split")
    @printf("  %-34s %8.1f %8.1f\n", "into air (r8)", mean(old_air), mean(new_air))
    @printf("  %-34s %8.1f %8.1f\n", "canopy PAR (r1)", mean(can_par), mean(can_par))
    @printf("  %-34s %8.1f %8.1f\n", "canopy NIR (r5)", mean(can_nir), mean(can_nir))
    @printf("  %-34s %8.1f %8.1f\n", "floor", 0.0, mean(flr))
    @printf("  %-34s %8.1f %8.1f\n", "cover (mostly lost outside)", 0.0, mean(cov))
    @printf("  %-34s %8.1f %8.1f   (as share of outside: %.2f -> %.2f)\n", "total inside (air+canopy+floor)",
            mean(old_air .+ can_par .+ can_nir), mean(new_air .+ can_par .+ can_nir .+ flr),
            mean(old_air .+ can_par .+ can_nir) / mean(I2), mean(new_air .+ can_par .+ can_nir .+ flr) / mean(I2))
end

okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g); okU = .!isnan.(U8obs); okU[1] = false
inset(ds) = in.(cd_, Ref(Set(ds)))
ALLD = inset(collect(0:114)); COM = inset(COMMON_TEST)
VP(T, RH) = RH ./ 100 .* Pws.(T) ./ 100; VPo = VP(Tobs_g, RHobs_g)
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
function hourly_rmse(v, o, mask)
    e = Float64[]
    for cols in values(HB); k = [c for c in cols if mask[c]]; length(k) >= 6 && push!(e, mean(v[k]) - mean(o[k])); end
    isempty(e) ? NaN : sqrt(mean(abs2, e))
end
PER = [(0, 8), (9, 48), (49, 100), (101, 109), (110, 114)]
famp(T5) = mean([maximum(T5[cd_ .== d]) - minimum(T5[cd_ .== d]) for d in 1:114])

runs = [("R0 D3 reference",        fl(pD3)),
        ("F1 floor_full",          fl(pD3; ff = 1.0)),
        ("K1 deep_k",              fl(pD3; dk = 1.0)),
        ("FK floor_full+deep_k",   fl(pD3; ff = 1.0, dk = 1.0)),
        ("S1 solar_split",         fl(pD3; ss = 1.0)),
        ("SF split+full+deep_k",   fl(pD3; ff = 1.0, ss = 1.0, dk = 1.0)),
        ("SFv SF, vent gain x0.5", update_params(fl(pD3; ff = 1.0, ss = 1.0, dk = 1.0), ["vc_gain"], [0.5 * p38.vc_gain]))]
res = Any[]
println("\nRunning $(length(runs)) full-season simulations ...")
for (nm_, p) in runs
    tt = @elapsed (us_ = simrun(p, weather_c, ctrl_v, u0_9))
    push!(res, (nm_, us_[2, :], rhf.(us_[2, :], us_[4, :]), clamp.(us_[7, :], 0.0, 1.0), us_[6, :], us_[9, :]))
    @printf("  %-24s %.1f s\n", nm_, tt)
end

println("\n", "="^168)
println("  FULL SEASON. hRMS = hourly T RMSE (K). Night bias (K) by period; day bias T (K), VP (hPa) and vent opening U8, all days.")
println("  T5 amp = mean daily floor-temperature range (K). Tdeep = deep-soil state (°C) on days 0/30/60/90/114.")
@printf("  %-24s | %6s %6s | %7s %7s %8s %7s %7s | %6s %6s %6s %7s %7s | %6s | %s\n", "run", "hRMS", "hRMS",
        "n 0-8", "n 9-48", "n 49-100", "n101-9", "n110-4", "T day", "VPday", "VPnit", "U8 day", "U8 nit", "T5amp", "Tdeep")
@printf("  %-24s | %6s %6s |\n", "", "all", "common")
println("  ", "-"^166)
for (nm_, T, RH, U, T5, Td) in res
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    tds = join([@sprintf("%5.1f", Td[col_of(d, 144)] - 273.15) for d in (0, 30, 60, 90, 114)], " ")
    @printf("  %-24s | %6.3f %6.3f | %+7.2f %+7.2f %+8.2f %+7.2f %+7.2f | %+6.2f %+6.2f %+6.2f %+7.3f %+7.3f | %6.2f | %s\n", nm_,
            hourly_rmse(T, Tobs_g, okT .& ALLD), hourly_rmse(T, Tobs_g, okT .& COM), nb...,
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& nite),
            mb(U, U8obs, okU .& ALLD .& dayt), mb(U, U8obs, okU .& ALLD .& nite), famp(T5), tds)
end

# ── k_ground scan ────────────────────────────────────────────────────────────
KG = [3.0, 5.0, 10.0, 20.0, p38.k_ground]
println("\n  k_ground scan: hourly T RMSE all days (night bias 49-100 / day bias) for each structure")
@printf("  %-24s |%s\n", "structure", join([@sprintf(" %21s", "k = $(round(k, digits = 1))") for k in KG], ""))
println("  ", "-"^(28 + 22 * length(KG)))
for (nm_, mk) in (("R0 D3 reference", kg -> fl(pD3; kg = kg)),
                  ("FK floor_full+deep_k", kg -> fl(pD3; ff = 1.0, dk = 1.0, kg = kg)),
                  ("SF split+full+deep_k", kg -> fl(pD3; ff = 1.0, ss = 1.0, dk = 1.0, kg = kg)))
    cells = String[]
    for kg in KG
        us_ = simrun(mk(kg), weather_c, ctrl_v, u0_9); T = us_[2, :]
        push!(cells, @sprintf(" %6.3f (%+5.2f/%+5.2f)", hourly_rmse(T, Tobs_g, okT .& ALLD),
                              mb(T, Tobs_g, okT .& nite .& inset(collect(49:100))), mb(T, Tobs_g, okT .& ALLD .& dayt)))
    end
    @printf("  %-24s |%s\n", nm_, join(cells, ""))
end
println("""
  How to read it:
   - R0 should give hRMS 1.254 (Note 44 D3); otherwise the flags-off path changed.
   - F1/FK: does closing the floor balance change the night bias pattern, and does the k_ground scan now have
     its best value BELOW the bound?
   - S1/SF: the solar partition removes ~1.0 x Iglob from the air. With the old parameters the day should turn
     cold and the vents should open less (U8 day bias falls). If VP day rises toward 0, the 'dry daytime air'
     (gamma4 at its bound) was caused by over-ventilating an over-heated air node.
   - SFv: if halving the vent gain recovers the day bias, the calibrated 2x vent gain (Notes 36-38) was
     compensating the extra solar heat.
   Parameters are the Note 38 posterior (old structure): this chooses the structure; the numbers come from the
   next calibration.""")

pl = plot(ylabel = "day bias (K)", xlabel = "day of season", title = "Daily DAY bias, model − measured", legend = :bottomleft,
          legendfontsize = 7, size = (1200, 450))
for (a, b) in TRAIN_BLOCKS; vspan!(pl, [a - 0.5, b + 0.5], color = :green, alpha = 0.08, label = ""); end
hline!(pl, [0], lc = :gray, ls = :dash, label = "")
for (k, (nm_, T, _, _, _, _)) in enumerate(res)
    k in (1, 4, 5, 6, 7) || continue
    plot!(pl, 0:114, [mb(T, Tobs_g, okT .& dayt .& (cd_ .== d)) for d in 0:114], lw = k == 1 ? 2.2 : 1.4, label = nm_)
end
pl2 = plot(ylabel = "night bias (K)", xlabel = "day of season", title = "Daily NIGHT bias, model − measured", legend = :bottomleft,
           legendfontsize = 7, size = (1200, 450))
for (a, b) in TRAIN_BLOCKS; vspan!(pl2, [a - 0.5, b + 0.5], color = :green, alpha = 0.08, label = ""); end
hline!(pl2, [0], lc = :gray, ls = :dash, label = "")
for (k, (nm_, T, _, _, _, _)) in enumerate(res)
    k in (1, 2, 4, 6, 7) || continue
    plot!(pl2, 0:114, [mb(T, Tobs_g, okT .& nite .& (cd_ .== d)) for d in 0:114], lw = k == 1 ? 2.2 : 1.4, label = nm_)
end
fig = joinpath(FIGDIR, "AiCU_floor_balance.png"); savefig(plot(pl, pl2, layout = (2, 1), size = (1200, 900)), fig)
println("Figure: ", fig)
