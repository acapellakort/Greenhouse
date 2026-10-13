# =============================================================================
# diagnose_heat_ctrl.jl -- test the default-off heating controller flag (heat_ctrl)
# =============================================================================
# Run from the GreenhouseSim folder (~10 s of simulation, no MCMC):
#     julia --project=. test/diagnose_heat_ctrl.jl
#
# 0. Bit-for-bit check: heat_ctrl = 0 passed explicitly must reproduce the base run exactly.
# Then, from N29 climate values with U8_min_leak = 0:
#   H0 base           heat_ctrl off, pipe input = measured / "setpoint + 25 K" estimate (as before)
#   H1 ctrl           heat_ctrl on, pipe input = measured or Tout + 2 (no estimate); Kp 10, tau_i 30 min
#   H2 ctrl dT -1     H1 with the setpoint 1 K lower (zone/sensor offset)
#   H3 ctrl Kp 20     H1 with stiffer control
#   H4 ctrl + est     heat_ctrl on, but keeping the old "setpoint + 25 K" pipe estimate as input
#   H5 ctrl+vent^0.5  H1 with effective roof opening U8^0.5 (does the controller remove the night
#                     cost that blocked the vent fix in Note 32-e?)
# Metrics: night/day bias by window, PRE-DAWN bias (04-08 h), regression slope on (T_in - T_out),
# and an approximate heating-demand indicator (mean pipe-air difference the controller asks for).
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = get(ENV, "N_FULL", "115")
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const HC_FIG = joinpath(FIGDIR, "AiCU_heat_ctrl")
hasproperty(params, :heat_ctrl) || error("constants_climate.json has no heat_ctrl entries")
@printf("\nJSON: heat_ctrl=%.0f hc_Kp=%.1f hc_Ki=%.2e hc_Tmax=%.1f C hc_dT=%.1f\n",
        params.heat_ctrl, params.hc_Kp, params.hc_Ki, params.hc_Tmax - 273.15, params.hc_dT)

# pipe input without the heating estimate: measured if valid, else Tout + 2 K
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
mkctrl(tp, u8) = ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(u8), zc, Z(U10v), Z(U11v), Z(U12v), Z(tp))
base_meas  = SimBase(params, weather, mkctrl(Tpipe_meas, U8v), t0)
base_mv05  = SimBase(params, weather, mkctrl(Tpipe_meas, sqrt.(U8v)), t0)

nV = QoI_V; θ0 = [4.46e-5, 44.0, 24.5, 0.0]
hcn = ["heat_ctrl", "hc_Kp", "hc_Ki", "hc_dT"]
hc(on, kp, taui, dT) = [on, kp, kp / taui, dT]

# ── 0. bit-for-bit check ─────────────────────────────────────────────────────
Ta = collect(forward_map(θ0, nV, base, days, tspan, teval; u0 = copy(u0))[2])
Tb = collect(forward_map(vcat(θ0, 0.0), vcat(nV, "heat_ctrl"), base, days, tspan, teval; u0 = copy(u0))[2])
@printf("\n[check 0] heat_ctrl = 0 vs V2.2 path: max |ΔT| = %.3e K  (%s)\n", maximum(abs.(Ta .- Tb)),
        Ta == Tb ? "bit-for-bit identical" : "NOT identical")

runs = [("H0 base",          nV,               θ0,                                 base),
        ("H1 ctrl",          vcat(nV, hcn),    vcat(θ0, hc(1.0, 10.0, 1800.0, 0.0)),  base_meas),
        ("H2 ctrl dT-1",     vcat(nV, hcn),    vcat(θ0, hc(1.0, 10.0, 1800.0, -1.0)), base_meas),
        ("H3 ctrl Kp20",     vcat(nV, hcn),    vcat(θ0, hc(1.0, 20.0, 1800.0, 0.0)),  base_meas),
        ("H4 ctrl+est",      vcat(nV, hcn),    vcat(θ0, hc(1.0, 10.0, 1800.0, 0.0)),  base),
        ("H5 ctrl+vent^.5",  vcat(nV, hcn),    vcat(θ0, hc(1.0, 10.0, 1800.0, 0.0)),  base_mv05)]
T2s = Vector{Vector{Float64}}(); I5s = Vector{Vector{Float64}}()
println("Running $(length(runs)) full-season forward simulations ...")
for (name, nm, x, b) in runs
    p = update_params(b.params, nm, x)
    t = @elapsed (us = GreenhouseSim._run(p, b, days, tspan, teval, copy(u0), GreenhouseSim.Tsit5()))
    size(us, 2) == NCOL || error("$name: solver stopped early")
    push!(T2s, us[2, :]); push!(I5s, us[5, :])
    @printf("  %-16s %.1f s\n", name, t)
end

coldays = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
hrs     = mod.((predtimes .- t0) ./ 3600, 24)
okT   = .!isnan.(Tobs_g); okT[1] = false
nite  = A(Ig_v).(predtimes) .< 1.0
dayt  = A(Ig_v).(predtimes) .> 50.0
pred  = (hrs .>= 4) .& (hrs .< 8)
win(a, b) = (coldays .>= a) .& (coldays .<= b)
Toutg = A(Tout_v).(predtimes); scr = A(U1v).(predtimes); U11g = A(U11v).(predtimes)
mb(T, m) = any(m) ? mean(T[m] .- Tobs_g[m]) : NaN
rmse(T, m) = sqrt(mean((T[m] .- Tobs_g[m]) .^ 2))
mreg = okT .& nite
X = hcat(ones(count(mreg)), Tobs_g[mreg] .- Toutg[mreg], scr[mreg])

println("\n", "="^128)
@printf("  %-16s %7s %7s | %7s %7s | %7s %7s | %7s %7s | %7s %7s | %7s %6s | %8s\n", "run", "RMSE", "RMSE",
        "night", "night", "predawn", "predawn", "day", "day", "day", "", "slope", "", "ctrl")
@printf("  %-16s %7s %7s | %7s %7s | %7s %7s | %7s %7s | %7s %7s | %7s %6s | %8s\n", "", "0-41", "42-114",
        "0-41", "42-114", "0-41", "42-114", "0-10", "11-41", "42-114", "", "b", "c", "u5>0 %")
println("  ", "-"^126)
for (k, r) in enumerate(runs)
    T = T2s[k]; y = T[mreg] .- Tobs_g[mreg]; β = X \ y
    @printf("  %-16s %7.2f %7.2f | %+7.2f %+7.2f | %+7.2f %+7.2f | %+7.2f %+7.2f | %+7.2f %7s | %+7.3f %+6.2f | %7.0f%%\n", r[1],
            rmse(T, okT .& win(0, 41)), rmse(T, okT .& win(42, N_DAYS - 1)),
            mb(T, okT .& nite .& win(0, 41)), mb(T, okT .& nite .& win(42, N_DAYS - 1)),
            mb(T, okT .& pred .& win(0, 41)), mb(T, okT .& pred .& win(42, N_DAYS - 1)),
            mb(T, okT .& dayt .& win(0, 10)), mb(T, okT .& dayt .& win(11, 41)), mb(T, okT .& dayt .& win(42, N_DAYS - 1)), "",
            β[2], β[3], 100 * mean(I5s[k][2:end] .> 1e-6))
end
println("""
  slope b / c: night residual regression on (T_in - T_out) and screen, all nights (the controller should flatten b).
  ctrl u5>0 %: share of time the controller's integral is active (heating demand). Pre-dawn = 04-08 h.""")

# figure: nights/pre-dawn on a warm and two cold days + daily night bias
SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "21,30,77,100"), ","))
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
pn = plot(xlabel = "day", ylabel = "night bias [K]", title = "Daily night bias", legend = :bottomleft, size = (1100, 350))
for (k, r) in enumerate(runs); plot!(pn, dd, [mb(T2s[k], okT .& nite .& (coldays .== d)) for d in dd], lw = 2, label = r[1]); end
hline!(pn, [0], lc = :black, label = "")
savefig(plot(plot(ps..., layout = (1, length(ps))), pn, layout = grid(2, 1, heights = [0.55, 0.45]),
             size = (420 * length(ps), 760)), HC_FIG * ".png")
println("Figure: ", HC_FIG, ".png")
