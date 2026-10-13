# =============================================================================
# diagnose_closed_loop.jl -- Note 34 flags: vent_ctrl, pipe_lag, black_scr (forward runs)
# =============================================================================
# Run from the GreenhouseSim folder (~20 s of simulation, no MCMC):
#     julia --project=. test/diagnose_closed_loop.jl
#
# Starts from the Note 33 `blocks` posterior mean (all Note 32 flags on) and adds ONE
# Note 34 flag at a time. Full season, days 0-114.
#   0. Checks: (a) ControlInputsX with the new flags OFF and the 6-state u0 must be
#      bit-for-bit identical to ControlInputs; (b) the 8-state u0 with flags off is reported
#      (not bit-for-bit: extra zero-derivative states change the adaptive step norm).
#   C0 open loop      Note 33 blocks posterior mean, logged vents replayed
#   C1 vent_ctrl 10m  closed-loop vents, tau 10 min (gain 0.12/K, start -2 K, hum 0.12)
#   C2 vent_ctrl 30m  same, tau 30 min
#   C3 pipe_lag 60m   heating pipe as a state, tau 60 min (logged vents)
#   C4 black_scr .5   blackout screen adds night insulation, k_black 0.5 (logged vents)
#   C5 all            vent_ctrl 10 min + pipe_lag 60 min + black_scr 0.5
# Metrics: T RMSE/bias by set (blocks design: TRAIN/VAL/TEST/COMMON TEST), night/day bias,
# spike overshoot (days 11-41), and -- for closed-loop runs -- the VENT OBSERVABLE:
# hourly-mean model opening vs logged opening (daytime RMSE, bias).
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = "115"
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const CL_FIG = joinpath(FIGDIR, "AiCU_closed_loop")
for k in (:vent_ctrl, :pipe_lag, :black_scr, :heat_ctrl, :vapour_fix, :cond_heat)
    hasproperty(params, k) || error("constants_climate.json has no $k entry")
end

# ── inputs as in run_climate_calibration.jl, plus the ControlInputsX channels ──
srcv = Int[]; VentSPv = Float64[]; HumSPv = Float64[]; Blackv = Float64[]; U8log = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    To = _num(mt.Tout[i]); To_val = isnan(To) ? 15.0 : To
    pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    hv = filter(!isnan, [_num(vip.HeatTemp_Vip_1[i]), _num(vip.HeatTemp_Vip_2[i]), _num(vip.HeatTemp_Vip_3[i]),
                         _num(vip.HeatTemp_Vip_5[i]), _num(vip.HeatTemp_Vip_6[i])])
    Th = isempty(hv) ? 20.0 : mean(hv)
    valid = (!isnan(pg) && pg > PIPE_MIN_VALID) || (!isnan(pl) && pl > PIPE_MIN_VALID)
    push!(srcv, valid ? 1 : (To_val < Th - 2.0 && To_val < PIPE_WARM_TOUT) ? 2 : 3)
    vs = filter(!isnan, vcat([_num(vip[i, Symbol("VentLeew_Vip_$z")]) for z in 1:6],
                             [_num(vip[i, Symbol("VentWind_Vip_$z")]) for z in 1:6]))
    push!(VentSPv, (isempty(vs) ? 25.0 : mean(vs)) + 273.15)
    hd = _num(vip.HumDef_Vip[i]); push!(HumSPv, isnan(hd) ? 0.0 : hd)
    bs = _num(gh.BlackScr[i]); push!(Blackv, isnan(bs) ? 0.0 : clamp(bs / 100, 0.0, 1.0))
    vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i])
    push!(U8log, (isnan(vl) || isnan(vw)) ? NaN : clamp((vl + vw) / 200, 0.0, 1.0))
end
length(srcv) == length(ts) || error("row filter mismatch")
Tpipe_meas = [srcv[i] == 2 ? Tout_v[i] + 2.0 : Tpipev[i] for i in eachindex(Tpipev)]
lai_d = vcat(0.0, Float64.(cm_days_all)); lai_v = vcat(0.5, cm_L_all)
lai_of_day(x) = x <= lai_d[1] ? lai_v[1] : x >= lai_d[end] ? lai_v[end] :
    (i = searchsortedlast(lai_d, x); lai_v[i] + (x - lai_d[i]) / (lai_d[i+1] - lai_d[i]) * (lai_v[i+1] - lai_v[i]))
LAI_row = [lai_of_day((t - t0) / 86400.0) for t in ts]
weather_c = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, _t -> 18.0 + 273.15,
                          A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), A(LAI_row))
U8z = [isnan(v) ? 0.0 : v for v in U8log]
ctrl_plain = ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas))
ctrl_x     = ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                            Z(VentSPv), Z(HumSPv), Z(Blackv))

# Note 33 blocks posterior mean (model parameters only)
cfile = joinpath(POSTDIR, "AiCU_climcal_blocks_seed1.jld2")
isfile(cfile) || error("missing $cfile (Note 33 blocks chain)")
cc = JLD2.load(cfile); CQm = ["gamma4", "vent_exp", "vent_scr_perm", "hc_dT", "f_cond", "k_ground"]
pmean = vec(mean(cc["chain"][cc["burn_in"]+1:end, 1:length(CQm)], dims = 1))
println("\nNote 33 blocks posterior mean: ", join(["$q=$(round(v, sigdigits=3))" for (q, v) in zip(CQm, pmean)], "  "))
FIXN = ["nu4", "U8_min_leak", "heat_ctrl", "hc_Kp", "hc_Ki", "vent_screen", "stomata_fix", "vapour_fix", "cond_heat"]
FIXV = [4.46e-5, 0.0, 1.0, 10.0, 10.0 / 1800, 1.0, 1.0, 1.0, 1.0]
p33 = update_params(update_params(params, FIXN, FIXV), CQm, pmean)

i0r = findfirst(!isnan, U8log)
u0_6 = copy(u0); u0_8 = vcat(u0, U8log[i0r], u0[2])
simrun(p, ctrl, u0x) = GreenhouseSim._run(p, SimBase(p, weather_c, ctrl, t0), days, tspan, teval, copy(u0x), GreenhouseSim.Tsit5())

# ── 0. checks ────────────────────────────────────────────────────────────────
us_a = simrun(p33, ctrl_plain, u0_6); us_b = simrun(p33, ctrl_x, u0_6)
@printf("\n[check 0a] ControlInputsX, flags off, 6 states vs ControlInputs: max |ΔT| = %.3e K (%s)\n",
        maximum(abs.(us_a[2, :] .- us_b[2, :])), us_a[2, :] == us_b[2, :] ? "bit-for-bit identical" : "NOT identical")
us_c = simrun(p33, ctrl_x, u0_8)
@printf("[check 0b] 8 states, flags off vs 6 states: max |ΔT| = %.3e K, mean |ΔT| = %.3e K (step-size effect only)\n",
        maximum(abs.(us_a[2, :] .- us_c[2, :])), mean(abs.(us_a[2, :] .- us_c[2, :])))

VC = ["vent_ctrl", "vc_tau"]; PL = ["pipe_lag", "pl_tau"]; BS = ["black_scr", "k_black"]
runs = [("C0 open loop",    p33, ctrl_plain, u0_6, false),
        ("C1 vent_ctrl 10m", update_params(p33, VC, [1.0, 600.0]),  ctrl_x, u0_8, true),
        ("C2 vent_ctrl 30m", update_params(p33, VC, [1.0, 1800.0]), ctrl_x, u0_8, true),
        ("C3 pipe_lag 60m",  update_params(p33, PL, [1.0, 3600.0]), ctrl_x, u0_8, false),
        ("C4 black_scr .5",  update_params(p33, BS, [1.0, 0.5]),    ctrl_x, u0_6, false),
        ("C5 all",           update_params(p33, vcat(VC, PL, BS), [1.0, 600.0, 1.0, 3600.0, 1.0, 0.5]), ctrl_x, u0_8, true)]
T2s = Vector{Vector{Float64}}(); U8s = Vector{Vector{Float64}}()
println("Running $(length(runs)) full-season simulations ...")
for (name, p, ct, u0x, closed) in runs
    t = @elapsed (us = simrun(p, ct, u0x))
    size(us, 2) == NCOL || error("$name: solver stopped early")
    push!(T2s, us[2, :]); push!(U8s, closed ? clamp.(us[7, :], 0.0, 1.0) : A(U8z).(predtimes))
    @printf("  %-18s %.1f s\n", name, t)
end

# ── metrics ──────────────────────────────────────────────────────────────────
wk(a, b) = collect(a:b)
TRAIN = vcat(wk(7, 13), wk(28, 34), wk(63, 69), wk(91, 97)); VAL = vcat(wk(42, 48), wk(77, 83))
TEST = setdiff(0:114, TRAIN, VAL); COMMON = vcat(wk(49, 62), wk(70, 76), wk(84, 90), wk(98, 114))
cd_ = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
okT = .!isnan.(Tobs_g); okT[1] = false
nite = A(Ig_v).(predtimes) .< 1.0; dayt = A(Ig_v).(predtimes) .> 50.0
U8obs = A(U8log).(predtimes)
inset(ds) = in.(cd_, Ref(Set(ds)))
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
rmse(v, o, m) = any(m) ? sqrt(mean((v[m] .- o[m]) .^ 2)) : NaN
# hourly means for the vent observable
hb = [c for c in 2:NCOL]; hour_id = [(cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)) for c in hb]
function hourly(v, mask)
    acc = Dict{Tuple{Int,Int},Vector{Float64}}()
    for (k, c) in enumerate(hb); mask[c] || continue; push!(get!(acc, hour_id[k], Float64[]), v[c]); end
    acc
end
function vent_obs(U8m, ds)
    m = dayt .& .!isnan.(U8obs) .& inset(ds)
    hm = hourly(U8m, m); ho = hourly(U8obs, m)
    ks = [k for k in keys(hm) if length(hm[k]) >= 6]
    isempty(ks) && return (NaN, NaN)
    dm = [mean(hm[k]) - mean(ho[k]) for k in ks]
    sqrt(mean(abs2, dm)), mean(dm)
end
function overshoot(T)
    ov = Float64[]
    for d in 11:41
        m = okT .& dayt .& (cd_ .== d); count(m) < 12 && continue
        push!(ov, maximum(T[m]) - maximum(Tobs_g[m]))
    end
    mean(ov), count(>(8), ov)
end

println("\n", "="^120)
@printf("  %-18s | %-27s | %-27s | %-11s | %-21s\n", "", "COMMON TEST: RMSE bias nite day", "TEST: RMSE bias nite day", "spikes 11-41", "VENT obs (day, hourly)")
@printf("  %-18s | %6s %6s %6s %6s | %6s %6s %6s %6s | %6s %4s | %6s %6s %7s\n", "run", "RMSE", "bias", "night", "day",
        "RMSE", "bias", "night", "day", "ovrsh", ">8K", "RMSE", "bias", "(COMMON)")
println("  ", "-"^118)
for (k, r) in enumerate(runs)
    T = T2s[k]; mc = okT .& inset(COMMON); mt_ = okT .& inset(TEST)
    ov, n8 = overshoot(T); vr, vb = vent_obs(U8s[k], COMMON)
    @printf("  %-18s | %6.2f %+6.2f %+6.2f %+6.2f | %6.2f %+6.2f %+6.2f %+6.2f | %+6.2f %4d | %6.3f %+6.3f\n", r[1],
            rmse(T, Tobs_g, mc), mb(T, Tobs_g, mc), mb(T, Tobs_g, mc .& nite), mb(T, Tobs_g, mc .& dayt),
            rmse(T, Tobs_g, mt_), mb(T, Tobs_g, mt_), mb(T, Tobs_g, mt_ .& nite), mb(T, Tobs_g, mt_ .& dayt),
            ov, n8, vr, vb)
end
println("""
  VENT obs: hourly-mean opening, model minus logged, daytime, COMMON TEST days (open-loop rows replay the log -> ~0).
  For closed-loop rows this is the new observable: a model that needs much more (bias > 0) or less vent than the
  grower used to hold the same setpoint has an energy-balance error the temperature alone no longer shows.""")

SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "24,30,41,77"), ","))
ps = Any[]
for d in SHOWD
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    pT = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 3, label = "measured", title = "day $d", ylabel = "T air (°C)",
              legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6)
    for (k, r) in enumerate(runs); k in (1, 2, 3, 6) && plot!(pT, h, T2s[k][c] .- 273.15, lw = 1.3, label = r[1]); end
    pV = plot(h, U8obs[c], lc = :black, lw = 2, label = "logged vent", xlabel = "hour", ylabel = "U8", ylims = (-0.05, 1.05),
              legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6)
    for (k, r) in enumerate(runs); k in (2, 3, 6) && plot!(pV, h, U8s[k][c], lw = 1.3, label = r[1]); end
    push!(ps, pT); push!(ps, pV)
end
order = [ps[2(i-1)+r] for r in 1:2 for i in 1:length(SHOWD)]
savefig(plot(order..., layout = (2, length(SHOWD)), size = (420 * length(SHOWD), 650)), CL_FIG * ".png")
println("Figure: ", CL_FIG, ".png")
