# =============================================================================
# run_climate_calibration.jl -- climate-only t-walk with the Note 32 structural flags ON
# =============================================================================
# Run from the GreenhouseSim folder:
#   DESIGN=blocks julia --project=. test/run_climate_calibration.jl smoke      # ALWAYS FIRST
#   DESIGN=early  julia --project=. test/run_climate_calibration.jl smoke
#   DESIGN=blocks SEED=1 CHAIN=100000 BURN=40000 julia --project=. test/run_climate_calibration.jl full
#
# Engine flags ON (all default-off in constants_climate.json, set here):
#   heat_ctrl (Kp 10, Ki 10/1800), vent_screen, stomata_fix, vapour_fix, cond_heat, LAI(t) input,
#   pipe input = measured or Tout+2 (no "setpoint + 25 K" estimate), U8_min_leak = 0, nu4 = 4.46e-5.
# Free (8): gamma4, vent_exp, vent_scr_perm, hc_dT, f_cond, k_ground, sigma_T, sigma_RH.
# Likelihood: hourly bin means of T_air and RH over the TRAIN days (Gaussian, sigmas inferred).
#
# DESIGN (Note 32-j; both are run and compared):
#   blocks : TRAIN = days 7-13, 28-34, 63-69, 91-97 (warm .. cold/heated), VAL = 42-48, 77-83,
#            TEST = all other days
#   early  : TRAIN = days 7-13, 21-27, 35-41 (forecast framing: nothing after day 41), VAL = 28-34,
#            TEST = days 42-114
#   COMMON TEST (reported for both): days 49-62, 70-76, 84-90, 98-114 -- TEST in both designs.
# Each TRAIN/VAL block is simulated separately with 1 spin-up day (not scored) before it,
# starting from the measured air state.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = "115"
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))   # season inputs (ts, U*v, Tair_obs, ...)

const CMODE  = isempty(ARGS) ? "smoke" : ARGS[1]
CMODE in ("smoke", "full") || error("usage: ... run_climate_calibration.jl [smoke|full]")
const DESIGN = get(ENV, "DESIGN", "blocks")
DESIGN in ("blocks", "early") || error("DESIGN must be blocks or early")
const CSEED  = parse(Int, get(ENV, "SEED", "1")); Random.seed!(CSEED)
const CCHAIN = parse(Int, get(ENV, "CHAIN", "100000"))
const CBURN  = parse(Int, get(ENV, "BURN", "40000"))
const CSMOKE = 200
const CTAG   = "AiCU_climcal_$(DESIGN)_seed$(CSEED)"

wk(a, b) = collect(a:b)
const TRAIN_BLOCKS = DESIGN == "blocks" ? [(7, 13), (28, 34), (63, 69), (91, 97)] : [(7, 13), (21, 27), (35, 41)]
const VAL_BLOCKS   = DESIGN == "blocks" ? [(42, 48), (77, 83)] : [(28, 34)]
const TRAIN_DAYS   = reduce(vcat, [wk(a, b) for (a, b) in TRAIN_BLOCKS])
const VAL_DAYS     = reduce(vcat, [wk(a, b) for (a, b) in VAL_BLOCKS])
const TEST_DAYS    = DESIGN == "blocks" ? setdiff(0:114, TRAIN_DAYS, VAL_DAYS) : collect(42:114)
const COMMON_TEST  = vcat(wk(49, 62), wk(70, 76), wk(84, 90), wk(98, 114))

# ── inputs with the structural options ──────────────────────────────────────
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
lai_d = vcat(0.0, Float64.(cm_days_all)); lai_v = vcat(0.5, cm_L_all)
lai_of_day(x) = x <= lai_d[1] ? lai_v[1] : x >= lai_d[end] ? lai_v[end] :
    (i = searchsortedlast(lai_d, x); lai_v[i] + (x - lai_d[i]) / (lai_d[i+1] - lai_d[i]) * (lai_v[i+1] - lai_v[i]))
LAI_row = [lai_of_day((t - t0) / 86400.0) for t in ts]          # precomputed step input (fast)
weather_c  = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, _t -> 18.0 + 273.15,
                           A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), A(LAI_row))
controls_c = ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8v), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas))
FIXN = ["nu4", "U8_min_leak", "heat_ctrl", "hc_Kp", "hc_Ki", "vent_screen", "stomata_fix", "vapour_fix", "cond_heat"]
FIXV = [4.46e-5, 0.0, 1.0, 10.0, 10.0 / 1800, 1.0, 1.0, 1.0, 1.0]
params_c = update_params(params, FIXN, FIXV)
const base_c = SimBase(params_c, weather_c, controls_c, t0)

# ── free parameters ──────────────────────────────────────────────────────────
const CQ_MODEL = ["gamma4", "vent_exp", "vent_scr_perm", "hc_dT", "f_cond", "k_ground"]
const CQ       = vcat(CQ_MODEL, ["sigma_T", "sigma_RH"])
const NM       = length(CQ_MODEL)
const CQ_dict = Dict(                     # [min, centre, max]
    "gamma4"        => [20.0, 82.0, 300.0],
    "vent_exp"      => [0.3,  0.7,  1.5],
    "vent_scr_perm" => [0.0,  0.2,  1.0],
    "hc_dT"         => [-3.0, 0.0,  3.0],
    "f_cond"        => [0.0,  0.5,  1.0],
    "k_ground"      => [1.0,  10.0, 30.0],
    "sigma_T"       => [0.2,  1.5,  6.0],
    "sigma_RH"      => [0.5,  6.0,  25.0])
CSupp(z) = all(CQ_dict[CQ[k]][1] < z[k] < CQ_dict[CQ[k]][3] for k in eachindex(z))

# ── blocks: 1 spin-up day + scored days, measured initial state ─────────────
struct Block
    days::Vector{Float64}; u0::Vector{Float64}
    binsT::Vector{Vector{Int}}; obsT::Vector{Float64}
    binsR::Vector{Vector{Int}}; obsR::Vector{Float64}
end
const AT = A(Tair_obs); const AR = A(RH_obs); const AC = A(CO2_obs)
function make_block(a, b)
    d0 = max(a - 1, 0); ds = collect(Float64, d0:b)
    tstart = t0 + 86400.0 * d0
    Ta = AT(tstart); Ta = isnan(Ta) ? Tair_obs[findfirst(!isnan, Tair_obs)] : Ta
    Rh = AR(tstart); Rh = isnan(Rh) ? 80.0 : Rh
    Co = AC(tstart); Co = isnan(Co) ? 575.0 : Co
    u0b = [Ta, Ta, Co, Rh / 100 * Pws(Ta), 0.0, Ta]
    times = vcat(tstart, [t0 + d * 86400.0 + tj for d in ds for tj in teval])
    To = AT.(times); Ro = AR.(times)
    bT = Vector{Int}[]; oT = Float64[]; bR = Vector{Int}[]; oR = Float64[]
    for (k, d) in enumerate(ds)
        Int(d) < a && continue                                  # spin-up day not scored
        for hb in 0:23
            cols = [2 + (k - 1) * NPD + j for j in 0:NPD-1 if min(floor(Int, j * 300 / 3600), 23) == hb]
            okT = [c for c in cols if !isnan(To[c])]; okR = [c for c in cols if !isnan(Ro[c])]
            length(okT) >= length(cols) / 2 && (push!(bT, okT); push!(oT, mean(To[okT])))
            length(okR) >= length(cols) / 2 && (push!(bR, okR); push!(oR, mean(Ro[okR])))
        end
    end
    return Block(ds, u0b, bT, oT, bR, oR)
end
const TB = [make_block(a, b) for (a, b) in TRAIN_BLOCKS]
const NBT = sum(length(b.obsT) for b in TB); const NBR = sum(length(b.obsR) for b in TB)

function block_sse(x, b::Block)
    out = try forward_map(x, CQ_MODEL, base_c, b.days, tspan, teval; u0 = copy(b.u0)) catch; return (Inf, Inf) end
    T2 = out[2]; RH = out[3]
    length(T2) == 1 + length(b.days) * NPD || return (Inf, Inf)
    sT = 0.0; sR = 0.0
    for k in eachindex(b.binsT); sT += (binmean(T2, b.binsT[k]) - b.obsT[k])^2; end
    for k in eachindex(b.binsR); sR += (binmean(RH, b.binsR[k]) - b.obsR[k])^2; end
    return (isfinite(sT) ? sT : Inf, isfinite(sR) ? sR : Inf)
end
function cenergy(z)
    x = z[1:NM]; σT = z[NM+1]; σR = z[NM+2]
    sT = 0.0; sR = 0.0
    for b in TB
        a, c = block_sse(x, b); (isfinite(a) && isfinite(c)) || return Inf
        sT += a; sR += c
    end
    return NBT * log(σT) + sT / (2σT^2) + NBR * log(σR) + sR / (2σR^2)
end

# progress + checkpoint (JTwalk calls Supp exactly once per iteration, +2 at setup)
mutable struct CProg; total::Int; calls::Int; t0::Float64; tl::Float64; tc::Float64; Umin::Float64; obj::Any; end
const CP = CProg(0, 0, 0.0, 0.0, 0.0, Inf, nothing)
const CPSEC = CMODE == "smoke" ? 10.0 : (stdout isa Base.TTY ? 30.0 : 300.0)
const CCKPT = joinpath(POSTDIR, CTAG * "_ckpt.jld2")
function creport(force = false)
    now = time(); (force || now - CP.tl >= CPSEC) || return; CP.tl = now
    it = clamp(CP.calls - 2, 0, CP.total); el = now - CP.t0; r = it > 0 ? el / it : NaN
    nb = round(Int, 30 * it / CP.total)
    line = @sprintf("[%s%s] %5.1f%%  %d/%d  elapsed %s  %.0f ms/it  ETA %s  Umin %.1f", "#"^nb, "-"^(30 - nb),
                    100it / CP.total, it, CP.total, fmt_t(el), 1e3r, fmt_t(r * (CP.total - it)), CP.Umin)
    stdout isa Base.TTY ? print("\r", line) : println(line); flush(stdout)
end
function ccheckpoint()
    (CMODE == "full" && CP.obj !== nothing) || return
    rows = min(clamp(CP.calls - 2, 0, CP.total), size(CP.obj.Output, 1)); rows < 2 && return
    JLD2.jldopen(CCKPT * ".tmp", "w") do f
        f["chain"] = CP.obj.Output[1:rows, :]; f["QoI"] = CQ; f["design"] = DESIGN; f["seed"] = CSEED
    end
    mv(CCKPT * ".tmp", CCKPT; force = true)
end
function csupp(z)
    CP.calls += 1; creport()
    if time() - CP.tc >= 1800; CP.tc = time(); ccheckpoint(); end
    return CSupp(z)
end
cenergy_t(z) = (U = cenergy(z); U < CP.Umin && (CP.Umin = U); U)

centre = [CQ_dict[q][2] for q in CQ]
function cjitter(c)
    for _ in 1:1000; x = c .* (1 .+ 0.05 .* randn(length(c))); x[4] = c[4] + 0.2 * randn(); CSupp(x) && return x; end
    error("no start point in support")
end

println("\n", "="^78)
@printf("Climate calibration -- design=%s  mode=%s  seed=%d\n", DESIGN, CMODE, CSEED)
println("  TRAIN days: ", TRAIN_BLOCKS, "   VAL: ", VAL_BLOCKS)
@printf("  TRAIN hourly bins: T %d, RH %d  (%d simulated days per step incl. spin-up)\n",
        NBT, NBR, sum(length(b.days) for b in TB))
println("="^78)

# ── evaluation on the full season ───────────────────────────────────────────
function full_eval(x)
    p = update_params(params_c, CQ_MODEL, x)
    us = GreenhouseSim._run(p, base_c, days, tspan, teval, copy(u0), GreenhouseSim.Tsit5())
    return us[2, :], rhf.(us[2, :], us[4, :])
end
function eval_table(T2, RH; label = "")
    cd_ = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
    okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g)
    nite = A(Ig_v).(predtimes) .< 1.0; dayt = A(Ig_v).(predtimes) .> 50.0
    VPo = RHobs_g ./ 100 .* Pws.(Tobs_g) ./ 100; VPm = RH ./ 100 .* Pws.(T2) ./ 100
    sets = [("TRAIN", TRAIN_DAYS), ("VAL", VAL_DAYS), ("TEST", TEST_DAYS), ("COMMON TEST", COMMON_TEST)]
    @printf("\n  %s\n  %-12s %6s %7s %7s %7s %7s %7s\n", label, "set", "RMSE", "bias", "night", "day", "VPday", "VPnite")
    for (nm_, ds) in sets
        m = okT .& in.(cd_, Ref(Set(ds))); mv_ = okV .& in.(cd_, Ref(Set(ds)))
        mb(v, o, k) = any(k) ? mean(v[k] .- o[k]) : NaN
        @printf("  %-12s %6.2f %+7.2f %+7.2f %+7.2f %+7.2f %+7.2f\n", nm_, sqrt(mean((T2[m] .- Tobs_g[m]) .^ 2)),
                mb(T2, Tobs_g, m), mb(T2, Tobs_g, m .& nite), mb(T2, Tobs_g, m .& dayt),
                mb(VPm, VPo, mv_ .& dayt), mb(VPm, VPo, mv_ .& nite))
    end
end

if CMODE == "smoke"
    cenergy(centre)
    te = @elapsed (U0 = cenergy(centre))
    @printf("\n[smoke 1] energy at prior centre: %.1f   one evaluation: %.0f ms\n", U0, 1e3te)
    isfinite(U0) || error("energy not finite at prior centre")
    obj = jtwalk(n = length(CQ), U = cenergy_t, Supp = csupp)
    CP.total = CSMOKE; CP.calls = 0; CP.t0 = CP.tl = CP.tc = time(); CP.obj = obj
    tw = @elapsed Run!(obj, T = CSMOKE, x0 = cjitter(centre), xp0 = cjitter(centre))
    creport(true); println()
    O = obj.Output; acc = mean(any(O[k+1, 1:end-1] .!= O[k, 1:end-1]) for k in 1:size(O, 1)-1)
    @printf("[smoke 2] %d steps in %.1f s (%.0f ms/step), moved-fraction %.2f, U %.1f -> %.1f\n",
            CSMOKE, tw, 1e3tw / CSMOKE, acc, O[1, end], O[end, end])
    @printf("          projected full run (%d steps): %.1f h\n", CCHAIN, tw / CSMOKE * CCHAIN / 3600)
    T2c, RHc = full_eval(centre[1:NM]); eval_table(T2c, RHc; label = "Full-season check at the PRIOR CENTRE (not calibrated):")
    println("\nSmoke done. Paste this output before launching `full`.")
    exit(0)
end

# ── full run ─────────────────────────────────────────────────────────────────
obj = jtwalk(n = length(CQ), U = cenergy_t, Supp = csupp)
x0c, xp0c = cjitter(centre), cjitter(centre)
CP.total = CCHAIN; CP.calls = 0; CP.t0 = CP.tl = CP.tc = time(); CP.obj = obj
@printf("Progress every %.0f s; checkpoint every 30 min -> %s\n", CPSEC, CCKPT)
@time Run!(obj, T = CCHAIN, x0 = x0c, xp0 = xp0c)
creport(true); println()
JLD2.@save joinpath(POSTDIR, CTAG * ".jld2") chain=obj.Output QoI=CQ burn_in=CBURN design=DESIGN seed=CSEED train=TRAIN_BLOCKS val=VAL_BLOCKS fixn=FIXN fixv=FIXV
println("Chain saved: ", joinpath(POSTDIR, CTAG * ".jld2"))

post = obj.Output[CBURN+1:end, 1:end-1]; pm = vec(mean(post, dims = 1))
println("\n--- posterior ---")
for (i, q) in enumerate(CQ)
    lo, _, hi = CQ_dict[q]; qs = quantile(post[:, i], [0.025, 0.975]); span = hi - lo
    flag = qs[1] - lo < 0.02span ? "AT LOWER BOUND" : hi - qs[2] < 0.02span ? "AT UPPER BOUND" : ""
    @printf("  %-14s mean %9.4g   95%% [%9.4g, %9.4g]  %s\n", q, pm[i], qs[1], qs[2], flag)
end
acc = mean(any(obj.Output[k+1, 1:end-1] .!= obj.Output[k, 1:end-1]) for k in 1:size(obj.Output, 1)-1)
@printf("moved-fraction %.3f\n", acc)
T2m, RHm = full_eval(pm[1:NM]); eval_table(T2m, RHm; label = "Full-season evaluation at the POSTERIOR MEAN:")
