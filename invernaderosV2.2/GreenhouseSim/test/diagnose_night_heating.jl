# =============================================================================
# diagnose_night_heating.jl -- heating, soil anchor and minimum vent as causes of the 3-phase night bias (Note 42)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1-2 min, no MCMC):
#     julia --project=. test/diagnose_night_heating.jl
#
# All runs at the Note 38 pooled posterior mean (PRIORS=n37, closed-loop vents), full season, reltol 1e-4.
#   R0  current model (Note 38)
#   H1  heat_ctrl = 0      : only the LOGGED pipe heats (logged 0 = heating off; no PI anywhere)
#   H2  hc_log_pipe = 1    : logged pipe where it is warm, PI only where no warm pipe is logged   [new flag]
#   H3  hc_dT = -1.5 K     : PI aims 1.5 K below the logged heating setpoint
#   S1  soil anchor = 7-day trailing mean of measured air T (input change, replaces the fixed 18 °C)
#   V1  vc_min_log = 1     : minimum vent opening = logged VentLeewMin setpoint / 200 (input)     [new flag]
#   C1  S1 + H2 + V1
#   C2-C4  C1 with the pipe length gamma1 (pipe -> air convection) x1.5, x2, x3
# Night bias is reported for the Note 41 periods: days 0-8, 9-48, 49-100, 101-109, 110-114.
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n37"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
for k in (:hc_log_pipe, :vc_min_log)
    hasproperty(params, k) || error("constants_climate.json has no $k entry (Note 42 flags missing)")
end

# ── Note 38 posterior mean ───────────────────────────────────────────────────
parts = Vector{Vector{Float64}}()
for s in (1, 2)
    d = JLD2.load(joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n37_seed$(s).jld2"))
    push!(parts, vec(mean(d["chain"][d["burn_in"]+1:end, 1:NM], dims = 1)))
end
pm = (parts[1] .+ parts[2]) ./ 2
p38 = update_params(params_cl, CQ_MODEL, pm)

# ── new inputs: minimum vent setpoint and a moving soil anchor ──────────────
VentMinv = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    vm = filter(!isnan, [_num(vip[i, Symbol("VentLeewMin_Vip_$z")]) for z in 1:6])
    push!(VentMinv, isempty(vm) ? 0.0 : clamp(mean(vm) / 200, 0.0, 1.0))
end
length(VentMinv) == length(ts) || error("row filter mismatch (VentMin)")
@printf("\nLogged minimum-vent input (U8 scale): mean %.3f, max %.3f, share of rows > 0: %.2f\n",
        mean(VentMinv), maximum(VentMinv), mean(VentMinv .> 0))

Tfill = copy(Tair_obs); for i in eachindex(Tfill); isnan(Tfill[i]) && (Tfill[i] = i > 1 ? Tfill[i-1] : 291.15); end
soil7 = similar(Tfill); j0 = 1; acc = 0.0
for i in eachindex(ts)
    global acc += Tfill[i]
    while ts[j0] < ts[i] - 7 * 86400; global acc -= Tfill[j0]; global j0 += 1; end
    soil7[i] = acc / (i - j0 + 1)
end
@printf("Soil anchor (7-day mean of measured air T): day 0 %.1f °C, day 50 %.1f °C, day 100 %.1f °C (fixed anchor: 18.0 °C)\n",
        soil7[1] - 273.15, A(soil7)(t0 + 50 * 86400.0) - 273.15, A(soil7)(t0 + 100 * 86400.0) - 273.15)
weather_s = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, A(soil7),
                          A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), A(LAI_row))
ctrl_v = ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                        Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv))

function simrun(p, w, ctrl, u0x)
    us = Matrix{Float64}(undef, length(u0x), 1 + length(days) * NPD); us[:, 1] = u0x
    u = copy(u0x); col = 2
    for d in days
        ctx = GreenhouseSim.SimContext(p, w, ctrl, t0 + 86400.0 * d, 86400.0 * d)
        sol = GreenhouseSim.solve(GreenhouseSim.ODEProblem(GreenhouseSim.rhs!, u, tspan, ctx), GreenhouseSim.Tsit5();
                                  saveat = teval, reltol = RELTOL)
        length(sol.t) == NPD || error("solver stopped on day $d")
        S = Array(sol); us[:, col:col+NPD-1] = S; col += NPD; u = S[:, end]
    end
    us
end
up(p, nn, vv) = update_params(p, nn, vv)

# check 0: VentMin input present but flag off == 16-argument ControlInputsX
ua = simrun(p38, weather_c, ctrl_x, u0_8); ub = simrun(p38, weather_c, ctrl_v, u0_8)
@printf("[check 0] ControlInputsX with VentMin (flag off) vs without: max |ΔT| = %.3e K (%s)\n",
        maximum(abs.(ua[2, :] .- ub[2, :])), ua[2, :] == ub[2, :] ? "bit-for-bit identical" : "NOT identical")

runs = [("R0 current (N38)",      p38,                                  weather_c, ctrl_x),
        ("H1 logged pipe only",   up(p38, ["heat_ctrl"], [0.0]),         weather_c, ctrl_x),
        ("H2 log pipe + PI gaps", up(p38, ["hc_log_pipe"], [1.0]),       weather_c, ctrl_x),
        ("H3 hc_dT = -1.5",       up(p38, ["hc_dT"], [-1.5]),            weather_c, ctrl_x),
        ("S1 soil = 7-day mean",  p38,                                   weather_s, ctrl_x),
        ("V1 logged min vent",    up(p38, ["vc_min_log"], [1.0]),        weather_c, ctrl_v),
        ("C1 S1+H2+V1",           up(p38, ["hc_log_pipe", "vc_min_log"], [1.0, 1.0]), weather_s, ctrl_v),
        ("C2 C1, pipe x1.5",      up(p38, ["hc_log_pipe", "vc_min_log", "gamma1"], [1.0, 1.0, 1.5 * p38.gamma1]), weather_s, ctrl_v),
        ("C3 C1, pipe x2",        up(p38, ["hc_log_pipe", "vc_min_log", "gamma1"], [1.0, 1.0, 2.0 * p38.gamma1]), weather_s, ctrl_v),
        ("C4 C1, pipe x3",        up(p38, ["hc_log_pipe", "vc_min_log", "gamma1"], [1.0, 1.0, 3.0 * p38.gamma1]), weather_s, ctrl_v)]
@printf("Pipe length gamma1 = %.2f m/m2 (convection scales with it; C2-C4 scale it)\n", p38.gamma1)
res = Any[]
println("Running $(length(runs)) full-season simulations ...")
for (nm_, p, w, ct) in runs
    tt = @elapsed (us_ = simrun(p, w, ct, u0_8))
    push!(res, (nm_, us_[2, :], rhf.(us_[2, :], us_[4, :]), clamp.(us_[7, :], 0.0, 1.0)))
    @printf("  %-22s %.1f s\n", nm_, tt)
end

# ── metrics ──────────────────────────────────────────────────────────────────
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

println("\n", "="^140)
println("  FULL SEASON. hRMS = hourly T RMSE (K). Night bias (K) by period. Day bias, VP night (hPa) and U8 night bias over all days.")
@printf("  %-22s | %6s %6s | %7s %7s %7s %7s %7s | %6s %7s %8s\n", "run", "hRMS", "hRMS", "n 0-8", "n 9-48", "n 49-100",
        "n101-9", "n110-4", "day", "VPnite", "U8 nite")
@printf("  %-22s | %6s %6s | %7s %7s %7s %7s %7s | %6s %7s %8s\n", "", "all", "common", "", "", "", "", "", "all", "all", "all")
println("  ", "-"^138)
for (nm_, T, RH, U) in res
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    @printf("  %-22s | %6.3f %6.3f | %+7.2f %+7.2f %+7.2f %+7.2f %+7.2f | %+6.2f %+7.2f %+8.3f\n", nm_,
            hourly_rmse(T, Tobs_g, okT .& ALLD), hourly_rmse(T, Tobs_g, okT .& COM), nb...,
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& nite), mb(U, U8obs, okU .& ALLD .& nite))
end
println("""
  How to read it:
   - H1/H2/H3 remove the warm middle (49-109) without new cold nights -> the PI heating adds heat the real system
     did not; take the matching flag into the next calibration.
   - S1 lifts the cold start (0-8) -> the fixed 18 °C soil anchor is the cause there (soil ladder or moving anchor).
   - V1 changes little -> the minimum-vent setpoint is real but small; keep it as an input anyway.
   - C1 shows whether the fixes add up. Parameters are still the Note 38 posterior (tuned with the old
     structure), so this is a sensitivity test, not a calibration.""")

# ── figure: daily night bias ─────────────────────────────────────────────────
pl = plot(ylabel = "night bias (K)", xlabel = "day of season", title = "Daily night bias, model − measured", legend = :bottomleft,
          legendfontsize = 7, size = (1200, 450))
for (a, b) in TRAIN_BLOCKS; vspan!(pl, [a - 0.5, b + 0.5], color = :green, alpha = 0.08, label = ""); end
hline!(pl, [0], lc = :gray, ls = :dash, label = "")
for (k, (nm_, T, _, _)) in enumerate(res)
    k in (1, 3, 5, 7, 9, 10) || continue
    nbd = [mb(T, Tobs_g, okT .& nite .& (cd_ .== d)) for d in 0:114]
    plot!(pl, 0:114, nbd, lw = k == 1 ? 2.2 : 1.4, label = nm_)
end
fig = joinpath(FIGDIR, "AiCU_night_heating.png"); savefig(pl, fig); println("Figure: ", fig)
