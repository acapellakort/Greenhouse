# =============================================================================
# diagnose_vent_flags.jl -- vent_exp and vent_screen flags, with the heating controller on
# =============================================================================
# Run from the GreenhouseSim folder (~10 s of simulation, no MCMC):
#     julia --project=. test/diagnose_vent_flags.jl
#
# 0. Bit-for-bit check: vent_exp = 1, vent_screen = 0 passed explicitly == V2.2 path.
# All runs: N29 climate values, U8_min_leak = 0, heat_ctrl ON (Kp 10, tau_i 30 min, dT 0),
# pipe input = measured or Tout + 2 (the "setpoint + 25 K" estimate is retired, Note 32-g).
#   V0 ctrl               reference (= H1 of diagnose_heat_ctrl.jl)
#   V1 +screen .2         vent_screen on, vent_scr_perm 0.2, linear vents
#   V2 +vent^.5           concave vents, no screen throttling (= H5)
#   V3 +vent^.5 +scr .2   both
#   V4 +vent^.7 +scr .2   milder concavity
#   V5 +vent^.5 +scr .05  tighter screen
# Metrics: RMSE, night / pre-dawn / day bias by window, spike-period overshoot (11-41),
# night regression slopes, heating-active share, vapour-pressure bias (day / night).
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = get(ENV, "N_FULL", "115")
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const VF_FIG = joinpath(FIGDIR, "AiCU_vent_flags")
for k in (:heat_ctrl, :vent_exp, :vent_screen, :vent_scr_perm)
    hasproperty(params, k) || error("constants_climate.json has no $k entry")
end

srcv = Int[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    To = _num(mt.Tout[i]); To_val = isnan(To) ? 15.0 : To
    pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    hv = filter(!isnan, [_num(vip.HeatTemp_Vip_1[i]), _num(vip.HeatTemp_Vip_2[i]), _num(vip.HeatTemp_Vip_3[i]),
                         _num(vip.HeatTemp_Vip_5[i]), _num(vip.HeatTemp_Vip_6[i])])
    Th = isempty(hv) ? 20.0 : mean(hv)
    valid = (!isnan(pg) && pg > PIPE_MIN_VALID) || (!isnan(pl) && pl > PIPE_MIN_VALID)
    push!(srcv, valid ? 1 : (To_val < Th - 2.0 && To_val < PIPE_WARM_TOUT) ? 2 : 3)
end
length(srcv) == length(ts) || error("row filter mismatch")
Tpipe_meas = [srcv[i] == 2 ? Tout_v[i] + 2.0 : Tpipev[i] for i in eachindex(Tpipev)]
base_meas = SimBase(params, weather,
                    ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8v), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas)), t0)

nV = QoI_V; θ0 = [4.46e-5, 44.0, 24.5, 0.0]
# ── 0. bit-for-bit check ─────────────────────────────────────────────────────
Ta = collect(forward_map(θ0, nV, base, days, tspan, teval; u0 = copy(u0))[2])
Tb = collect(forward_map(vcat(θ0, 1.0, 0.0), vcat(nV, "vent_exp", "vent_screen"), base, days, tspan, teval; u0 = copy(u0))[2])
@printf("\n[check 0] vent_exp=1, vent_screen=0 vs V2.2 path: max |ΔT| = %.3e K  (%s)\n", maximum(abs.(Ta .- Tb)),
        Ta == Tb ? "bit-for-bit identical" : "NOT identical")

nm = vcat(nV, "heat_ctrl", "hc_Kp", "hc_Ki", "hc_dT", "vent_exp", "vent_screen", "vent_scr_perm")
cfg(ve, vs, perm) = vcat(θ0, 1.0, 10.0, 10.0 / 1800, 0.0, ve, vs, perm)
runs = [("V0 ctrl",            cfg(1.0, 0.0, 0.2)),
        ("V1 +scr.2",          cfg(1.0, 1.0, 0.2)),
        ("V2 +vent^.5",        cfg(0.5, 0.0, 0.2)),
        ("V3 +vent^.5 +scr.2", cfg(0.5, 1.0, 0.2)),
        ("V4 +vent^.7 +scr.2", cfg(0.7, 1.0, 0.2)),
        ("V5 +vent^.5 +scr.05",cfg(0.5, 1.0, 0.05))]
T2s = Vector{Vector{Float64}}(); RHs = Vector{Vector{Float64}}(); I5s = Vector{Vector{Float64}}()
println("Running $(length(runs)) full-season forward simulations ...")
for (name, x) in runs
    p = update_params(params, nm, x)
    t = @elapsed (us = GreenhouseSim._run(p, base_meas, days, tspan, teval, copy(u0), GreenhouseSim.Tsit5()))
    size(us, 2) == NCOL || error("$name: solver stopped early")
    push!(T2s, us[2, :]); push!(RHs, rhf.(us[2, :], us[4, :])); push!(I5s, us[5, :])
    @printf("  %-20s %.1f s\n", name, t)
end

coldays = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
hrs   = mod.((predtimes .- t0) ./ 3600, 24)
okT   = .!isnan.(Tobs_g); okT[1] = false
okV   = okT .& .!isnan.(RHobs_g)
nite  = A(Ig_v).(predtimes) .< 1.0
dayt  = A(Ig_v).(predtimes) .> 50.0
pred  = (hrs .>= 4) .& (hrs .< 8)
win(a, b) = (coldays .>= a) .& (coldays .<= b)
Toutg = A(Tout_v).(predtimes); scr = A(U1v).(predtimes); U11g = A(U11v).(predtimes)
VPobs = RHobs_g ./ 100 .* Pws.(Tobs_g) ./ 100
mb(T, o, m) = any(m) ? mean(T[m] .- o[m]) : NaN
rmse(T, m) = sqrt(mean((T[m] .- Tobs_g[m]) .^ 2))
mreg = okT .& nite
X = hcat(ones(count(mreg)), Tobs_g[mreg] .- Toutg[mreg], scr[mreg])
function overshoot(T)
    ov = Float64[]
    for d in 11:41
        m = okT .& dayt .& (coldays .== d); count(m) < 12 && continue
        push!(ov, maximum(T[m]) - maximum(Tobs_g[m]))
    end
    mean(ov), count(>(8), ov)
end

println("\n", "="^146)
@printf("  %-20s %6s %6s | %6s %6s | %6s %6s | %6s %6s %6s | %6s %4s | %6s %6s | %5s | %6s %6s\n", "run",
        "RMSE", "RMSE", "night", "night", "prdwn", "prdwn", "day", "day", "day", "ovrsh", ">8K", "b", "c", "heat", "VPday", "VPnite")
@printf("  %-20s %6s %6s | %6s %6s | %6s %6s | %6s %6s %6s | %6s %4s | %6s %6s | %5s | %6s %6s\n", "",
        "0-41", "42+", "0-41", "42+", "0-41", "42+", "0-10", "11-41", "42+", "11-41", "", "", "", "%", "hPa", "hPa")
println("  ", "-"^144)
for (k, (name, _)) in enumerate(runs)
    T = T2s[k]; VP = RHs[k] ./ 100 .* Pws.(T) ./ 100
    β = X \ (T[mreg] .- Tobs_g[mreg]); ov, n8 = overshoot(T)
    @printf("  %-20s %6.2f %6.2f | %+6.2f %+6.2f | %+6.2f %+6.2f | %+6.2f %+6.2f %+6.2f | %+6.2f %4d | %+6.3f %+6.2f | %5.0f | %+6.2f %+6.2f\n",
            name, rmse(T, okT .& win(0, 41)), rmse(T, okT .& win(42, N_DAYS - 1)),
            mb(T, Tobs_g, okT .& nite .& win(0, 41)), mb(T, Tobs_g, okT .& nite .& win(42, N_DAYS - 1)),
            mb(T, Tobs_g, okT .& pred .& win(0, 41)), mb(T, Tobs_g, okT .& pred .& win(42, N_DAYS - 1)),
            mb(T, Tobs_g, okT .& dayt .& win(0, 10)), mb(T, Tobs_g, okT .& dayt .& win(11, 41)),
            mb(T, Tobs_g, okT .& dayt .& win(42, N_DAYS - 1)), ov, n8, β[2], β[3],
            100 * mean(I5s[k][2:end] .> 1e-6), mb(VP, VPobs, okV .& dayt), mb(VP, VPobs, okV .& nite))
end
println("""
  heat %: share of time the heating controller's integral is active. VP: vapour-pressure bias, all days.
  Target: night and pre-dawn ~0 (kept from V0), day bias 11-41 / 42+ toward 0, fewer >8 K days, heat % not
  pinned near 100 (that would mean the pipe is saturated and the night bias is being hidden).""")

SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "24,30,41,100"), ","))
ps = Any[]
for d in SHOWD
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    p = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 3, label = "measured", title = "day $d", xlabel = "hour",
             ylabel = "T air (°C)", legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6)
    plot!(p, h, U11g[c] .- 273.15, lc = :red, ls = :dot, label = "heating SP")
    for (k, r) in enumerate(runs); plot!(p, h, T2s[k][c] .- 273.15, lw = 1.3, label = r[1]); end
    push!(ps, p)
end
dd = 0:N_DAYS-1
pn = plot(xlabel = "day", ylabel = "bias [K]", title = "Daily NIGHT bias (solid) and DAY bias (dashed)",
          legend = :bottomleft, size = (1100, 350))
for (k, r) in enumerate(runs)
    plot!(pn, dd, [mb(T2s[k], Tobs_g, okT .& nite .& (coldays .== d)) for d in dd], lw = 2, label = r[1], c = k)
    plot!(pn, dd, [mb(T2s[k], Tobs_g, okT .& dayt .& (coldays .== d)) for d in dd], lw = 1, ls = :dash, label = "", c = k)
end
hline!(pn, [0], lc = :black, label = "")
savefig(plot(plot(ps..., layout = (1, length(ps))), pn, layout = grid(2, 1, heights = [0.55, 0.45]),
             size = (420 * length(ps), 760)), VF_FIG * ".png")
println("Figure: ", VF_FIG, ".png")
