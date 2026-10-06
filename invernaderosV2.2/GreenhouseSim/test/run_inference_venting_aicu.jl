# =============================================================================
# run_inference_venting_aicu.jl
#   t-walk MCMC for nu4 (leakage) + gamma4 (stomata) + k_ground (floor)
#   on AiCU AGC-2018 data with a disjoint TRAIN / VAL split.
#
# TRAIN window : days 21–23  (serial 43347–43349)  U8 = 0.670 / 0.503 / 0.317
#   Day 21 is the best-ventilated day in the series (U8=0.670) — maximum
#   identifiability of nu4.  Days 22–23 add weather diversity.
#
# VAL   window : days  6–8   (serial 43332–43334)  U8 = 0.582 / 0.561 / 0.551
#   Earlier week, different synoptic pattern. Held out during calibration.
#   VAL RMSE_T is the honest diagnostic number to compare against v3's 4.09 K.
#
# CAVEAT: both windows are within the same 42-day block (Aug–Sep 2018).
# This is within-dataset TRAIN/VAL, not the disjoint-year standard.
# Multi-year validation must wait for additional AiCU seasons.
#
# Tpipe: uses v3 heating-estimate (Tpipe = T_heat_setpoint + PIPE_DELTA_K when
# outdoor < setpoint − 2°C), same as rmse_vanthoor_diagnostic_v3.jl.
#
# Run from the GreenhouseSim folder:
#   julia --project=. test/run_inference_venting_aicu.jl
# Expected wall time: ~35–50 min  (120 K MCMC iters × 3-day forward map)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, JTwalk, Plots, Printf, Statistics

# ── paths ────────────────────────────────────────────────────────────────────
const INV    = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TEAM   = "AiCU"
const GH_CSV  = joinpath(INV, "invernaderos2.0", "data", TEAM, "Greenhouse_climate.csv")
const MET_CSV = joinpath(INV, "invernaderos2.0", "data", "meteo.csv")
const VIP_CSV = joinpath(INV, "invernaderos2.0", "data", TEAM, "vip.csv")
const CLM_JS  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
const CRP_JS  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
const FIGDIR  = joinpath(@__DIR__, "figures"); mkpath(FIGDIR)

# ── TRAIN / VAL windows ───────────────────────────────────────────────────────
const TRAIN_START = 43347.0   # AiCU day 21  (U8=0.670, best for nu4 identification)
const TRAIN_DAYS  = 3         # days 21-23
const VAL_START   = 43332.0   # AiCU day  6  (U8=0.582)
const VAL_DAYS    = 3         # days  6-8

# ── settings ──────────────────────────────────────────────────────────────────
const LAI_FIX      = 2.0
const DAY_LO       = 7.0     # scored daytime window [h]
const DAY_HI       = 20.0
const SIG_T        = 1.5     # observation noise [K]
const SIG_RH       = 8.0     # observation noise [%RH]
const CHAIN_LENGTH = 120_000
const BURN_IN      = 40_000
const PIPE_DELTA_K = 20.0    # v3: pipe above setpoint when heating [K]
const PIPE_MIN_VALID = 5.0   # pipe readings ≤ this = "off" [°C]

# ── parameters to calibrate ──────────────────────────────────────────────────
const QoI = ["nu4", "gamma4", "k_ground"]
const QoI_dict = Dict(
    "nu4"      => [1.0e-5, 1.0e-4, 1.0e-1],  # heat leakage
    "gamma4"   => [20.0,   82.0,   300.0],    # min stomatal resistance
    "k_ground" => [0.5,    5.0,    50.0],     # air↔floor conductance
)

# ── helpers ───────────────────────────────────────────────────────────────────
const EXCEL_UNIX0 = 25569.0
const SIGMA_SB    = 5.670374419e-8
serial2unix(s) = (s - EXCEL_UNIX0) * 86400.0
_num(x) = (x isa Number && !(x isa Missing) && !isnan(Float64(x))) ? Float64(x) : NaN
stepfun(ts, vs) = t -> @inbounds vs[clamp(searchsortedlast(ts, t), 1, length(ts))]
_q2(T) = T > 0 ? 611.21*exp((18.678-T/234.5)*(T/(257.14+T))) :
                  611.21*exp((23.036-T/333.7)*(T/(279.82+T)))
_vpout(Tc, RH) = _q2(Tc) * clamp(RH, 1.0, 100.0) / 100.0

# ── data loader (shared for TRAIN and VAL) ───────────────────────────────────
gh  = DataFrame(CSV.File(GH_CSV))
mt  = DataFrame(CSV.File(MET_CSV))
vip = DataFrame(CSV.File(VIP_CSV))
n   = min(nrow(gh), nrow(mt), nrow(vip))

function load_window(win_lo, ndays)
    win_hi = win_lo + ndays
    ts=Float64[]; Tout_v=Float64[]; Ig_v=Float64[]; Wind_v=Float64[]
    VPout_v=Float64[]; Tsky_v=Float64[]
    U8_v=Float64[]; U1_v=Float64[]; U12_v=Float64[]; Tpipe_v=Float64[]
    Tair_obs=Float64[]; RH_obs=Float64[]; CO2_obs=Float64[]

    for i in 1:n
        g = _num(gh.GHtime[i])
        (isnan(g) || g < win_lo || g >= win_hi) && continue
        Ta = _num(gh.Tair[i]); isnan(Ta) && continue

        To  = _num(mt.Tout[i]);  Igv = _num(mt.Iglob[i])
        Ws  = _num(mt.Windsp[i]); Rho = _num(mt.Rhout[i]); Py = _num(mt.Pyrgeo[i])
        vl  = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i])
        pg  = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
        en  = _num(gh.EnScr[i]); al = _num(gh.AssimLight[i])
        Co  = _num(gh.CO2air[i]); Rh = _num(gh.RHair[i])

        # zone-averaged heating setpoint from vip.csv
        h1 = _num(vip.HeatTemp_Vip_1[i]); h2 = _num(vip.HeatTemp_Vip_2[i])
        h3 = _num(vip.HeatTemp_Vip_3[i]); h5 = _num(vip.HeatTemp_Vip_5[i])
        h6 = _num(vip.HeatTemp_Vip_6[i])
        hvals = filter(!isnan, [h1,h2,h3,h5,h6])
        T_heat_C = isempty(hvals) ? 20.0 : mean(hvals)

        ToK  = (isnan(To) ? 15.0 : To) + 273.15
        Ldn  = SIGMA_SB * ToK^4 + (isnan(Py) ? 0.0 : Py)
        Tsk  = isnan(Py) ? (-0.4+273.15) : clamp((Ldn/SIGMA_SB)^0.25, ToK-30.0, ToK+2.0)

        # v3 Tpipe: measured if valid; else heating-setpoint estimate
        pg_c = (!isnan(pg) && pg > PIPE_MIN_VALID) ? pg : NaN
        pl_c = (!isnan(pl) && pl > PIPE_MIN_VALID) ? pl : NaN
        pipe_max = max(isnan(pg_c) ? -Inf : pg_c, isnan(pl_c) ? -Inf : pl_c)
        To_val = isnan(To) ? 15.0 : To
        if pipe_max == -Inf
            Tpipe = (To_val < T_heat_C - 2.0 ?
                     (T_heat_C + PIPE_DELTA_K) : (To_val + 2.0)) + 273.15
        else
            Tpipe = pipe_max + 273.15
        end

        push!(ts,      serial2unix(g))
        push!(Tout_v,  ToK)
        push!(Ig_v,    isnan(Igv) ? 0.0 : max(Igv, 0.0))
        push!(Wind_v,  isnan(Ws)  ? 2.0 : Ws)
        push!(VPout_v, _vpout(To_val, isnan(Rho) ? 80.0 : Rho))
        push!(Tsky_v,  Tsk)
        push!(U8_v,    clamp(((isnan(vl) ? 0.0 : vl)+(isnan(vw) ? 0.0 : vw))/2/100, 0.0, 1.0))
        push!(U1_v,    isnan(en) ? 0.0 : clamp(en/100, 0.0, 1.0))
        push!(U12_v,   isnan(al) ? 0.0 : clamp(al/100, 0.0, 1.0))
        push!(Tpipe_v, Tpipe)
        push!(Tair_obs, Ta + 273.15)
        push!(RH_obs,   isnan(Rh) ? NaN : Rh)
        push!(CO2_obs,  isnan(Co) ? 575.0 : Co*536.4/(Ta+273.15))
    end
    @assert length(ts) > 20 "No rows found for window lo=$win_lo"
    if median(Wind_v) > 25; Wind_v ./= 3.6; end
    return (ts=ts, Tout=Tout_v, Ig=Ig_v, Wind=Wind_v, VPout=VPout_v, Tsky=Tsky_v,
            U8=U8_v, U1=U1_v, U12=U12_v, Tpipe=Tpipe_v,
            Tair=Tair_obs, RH=RH_obs, CO2=CO2_obs)
end

function make_base(d)
    t0 = d.ts[1]
    A(v) = stepfun(d.ts, v)
    Z(v) = stepfun(d.ts .- t0, v)
    zc   = Z(zeros(length(d.ts)))
    params  = load_params(CLM_JS, CRP_JS)
    weather = WeatherInputs(A(d.Ig), A(d.Tsky), A(d.Tout),
                            _t->273.15, _t->18.0+273.15,
                            A(d.Wind), _t->834.7, A(d.VPout), A(d.Ig), _t->LAI_FIX)
    controls = ControlInputs(Z(d.U1), zc, zc, zc, zc, zc,
                             zc, Z(d.U8), zc, zc, zc, Z(d.U12), Z(d.Tpipe))
    SimBase(params, weather, controls, t0)
end

function prediction_grid(d, ndays)
    t0 = d.ts[1]
    days  = collect(0.0:(ndays-1))
    teval = 0.0:900.0:86400.0
    predtimes = vcat(t0, [t0 + day*86400.0 + tj for day in days for tj in teval])
    hours = mod.((predtimes .- t0)./3600.0, 24.0)
    A(v) = stepfun(d.ts, v)
    Tmeas  = A(d.Tair).(predtimes)
    RHmeas = A(d.RH).(predtimes)
    mask = (hours .>= DAY_LO) .& (hours .<= DAY_HI) .& .!isnan.(Tmeas) .& .!isnan.(RHmeas)
    (days=days, tspan=(0.0,86400.0), teval=teval,
     predtimes=predtimes, hours=hours,
     Tmeas=Tmeas, RHmeas=RHmeas, mask=mask)
end

# ── load windows ──────────────────────────────────────────────────────────────
println("Loading TRAIN window (days 21–23, serial $(TRAIN_START)–$(TRAIN_START+TRAIN_DAYS))…")
dTR = load_window(TRAIN_START, TRAIN_DAYS)
bTR = make_base(dTR)
gTR = prediction_grid(dTR, TRAIN_DAYS)
@printf("TRAIN scored points (daytime): %d of %d\n", sum(gTR.mask), length(gTR.mask))
@printf("TRAIN U8 mean=%.3f  Tpipe mean=%.1f°C\n",
    mean(dTR.U8), mean(dTR.Tpipe)-273.15)

println("\nLoading VAL window (days 6–8, serial $(VAL_START)–$(VAL_START+VAL_DAYS))…")
dVAL = load_window(VAL_START, VAL_DAYS)
bVAL = make_base(dVAL)
gVAL = prediction_grid(dVAL, VAL_DAYS)
@printf("VAL   scored points (daytime): %d of %d\n", sum(gVAL.mask), length(gVAL.mask))
@printf("VAL   U8 mean=%.3f  Tpipe mean=%.1f°C\n",
    mean(dVAL.U8), mean(dVAL.Tpipe)-273.15)

# ── initial state ─────────────────────────────────────────────────────────────
function u0_from(d)
    V1_0 = (isnan(d.RH[1]) ? 80.0 : d.RH[1]) / 100 * Pws(d.Tair[1])
    [d.Tair[1], d.Tair[1], d.CO2[1], V1_0, 0.0, 18.0+273.15]
end
u0_TR  = u0_from(dTR)
u0_VAL = u0_from(dVAL)

# ── MCMC on TRAIN ─────────────────────────────────────────────────────────────
PriorSupp(x) = all(QoI_dict[QoI[k]][1] < x[k] < QoI_dict[QoI[k]][3] for k in eachindex(x))

function energy_train(x)
    _, T2, RH, _ = forward_map(x, QoI, bTR, gTR.days, gTR.tspan, gTR.teval; u0=copy(u0_TR))
    Tp = collect(T2); RHp = collect(RH)
    logL = -sum(((gTR.Tmeas[gTR.mask] .- Tp[gTR.mask])./SIG_T).^2) -
            sum(((gTR.RHmeas[gTR.mask] .- RHp[gTR.mask])./SIG_RH).^2)
    return -logL
end

n_p = length(QoI)
obj = jtwalk(n=n_p, U=energy_train, Supp=PriorSupp)
rng_init() = [rand()*(QoI_dict[q][3]-QoI_dict[q][1])+QoI_dict[q][1] for q in QoI]
x0  = rng_init(); xp0 = rng_init()

println("\nInferring: ", QoI)
println("Chain: $CHAIN_LENGTH iters, burn-in: $BURN_IN")
println("Estimated wall time: ~35–50 min")
flush(stdout)
@time Run!(obj, T=CHAIN_LENGTH, x0=x0, xp0=xp0)

# ── posteriors ────────────────────────────────────────────────────────────────
println("\n─── TRAIN posteriors (post burn-in) ──────────────────────────────")
med = Float64[]
for (i, name) in enumerate(QoI)
    c = obj.Output[BURN_IN:end, i]
    q = quantile(c, [0.025, 0.5, 0.975]); push!(med, q[2])
    @printf("%-10s  median=%.3e   95%% CI [%.3e, %.3e]   prior-center=%.3e\n",
            name, q[2], q[1], q[3], QoI_dict[name][2])
    p1 = plot(obj.Output[:,i], title="chain: $name", label="")
    vline!(p1, [BURN_IN], lc=:red, ls=:dash, label="burn-in")
    p2 = histogram(c, title="posterior: $name", label="")
    vline!(p2, [QoI_dict[name][2]], lc=:green, lw=3, label="prior centre")
    savefig(p1, joinpath(FIGDIR, "aicu_inf_$(name)_chain.png"))
    savefig(p2, joinpath(FIGDIR, "aicu_inf_$(name)_posterior.png"))
end

pE = plot(obj.Output[:,end], title="Energy (−logL)", xlabel="iter", label="")
vline!(pE, [BURN_IN], lc=:red, ls=:dash, label="burn-in")
savefig(pE, joinpath(FIGDIR, "aicu_inf_energy.png"))

# ── TRAIN fit ────────────────────────────────────────────────────────────────
_, T2_tr, RH_tr, _ = forward_map(med, QoI, bTR, gTR.days, gTR.tspan, gTR.teval; u0=copy(u0_TR))
ph_tr = (gTR.predtimes .- dTR.ts[1])./3600.0
Tp_tr = collect(T2_tr) .- 273.15; Tm_tr = gTR.Tmeas .- 273.15
rmse_tr = sqrt(mean((Tp_tr[gTR.mask] .- Tm_tr[gTR.mask]).^2))
bias_tr = mean(Tp_tr[gTR.mask] .- Tm_tr[gTR.mask])

pT_tr = plot(ph_tr, Tm_tr, lc=:black, lw=2, label="measured", xlabel="hours", ylabel="T air (°C)")
plot!(pT_tr, ph_tr, Tp_tr, lc=:orange, lw=2, label="model (posterior median)")
vspan!(pT_tr, [DAY_LO, DAY_HI], alpha=0.08, lc=:blue, label="scored window")
title!(pT_tr, "TRAIN fit  RMSE=$(round(rmse_tr,digits=2)) K  bias=$(round(bias_tr,digits=2)) K")
savefig(pT_tr, joinpath(FIGDIR, "aicu_inf_train_Tfit.png"))

# ── VAL evaluation ────────────────────────────────────────────────────────────
println("\n─── VAL evaluation (held-out, days 6–8) ─────────────────────────")
_, T2_val, RH_val, _ = forward_map(med, QoI, bVAL, gVAL.days, gVAL.tspan, gVAL.teval; u0=copy(u0_VAL))
ph_val = (gVAL.predtimes .- dVAL.ts[1])./3600.0
Tp_val = collect(T2_val) .- 273.15; Tm_val = gVAL.Tmeas .- 273.15
rmse_val = sqrt(mean((Tp_val[gVAL.mask] .- Tm_val[gVAL.mask]).^2))
bias_val = mean(Tp_val[gVAL.mask] .- Tm_val[gVAL.mask])
@printf("VAL  RMSE_T = %.2f K   bias = %+.2f K\n", rmse_val, bias_val)

pT_val = plot(ph_val, Tm_val, lc=:black, lw=2, label="measured", xlabel="hours", ylabel="T air (°C)")
plot!(pT_val, ph_val, Tp_val, lc=:steelblue, lw=2, label="model (posterior median)")
vspan!(pT_val, [DAY_LO, DAY_HI], alpha=0.08, lc=:blue, label="scored window")
title!(pT_val, "VAL fit  RMSE=$(round(rmse_val,digits=2)) K  bias=$(round(bias_val,digits=2)) K")
savefig(pT_val, joinpath(FIGDIR, "aicu_inf_val_Tfit.png"))

# ── summary ───────────────────────────────────────────────────────────────────
println("\n═══════════════════════════════════════════════════════════════════")
println("  Inference summary — AiCU, v3 Tpipe, TRAIN days 21-23")
println("═══════════════════════════════════════════════════════════════════")
for (i,name) in enumerate(QoI)
    c = obj.Output[BURN_IN:end,i]
    q = quantile(c, [0.025,0.5,0.975])
    @printf("  %-10s  %.3e  [%.3e, %.3e]\n", name, q[2], q[1], q[3])
end
@printf("\n  TRAIN RMSE_T : %.2f K  bias : %+.2f K  (days 21-23, scored 7-20h)\n", rmse_tr, bias_tr)
@printf("  VAL   RMSE_T : %.2f K  bias : %+.2f K  (days  6-8, scored 7-20h)\n", rmse_val, bias_val)
println("\n  Note: within-dataset TRAIN/VAL only (42-day block Aug–Sep 2018).")
println("  Multi-year validation pending additional AiCU seasons.")
println("═══════════════════════════════════════════════════════════════════")
println("\nFigures in: ", FIGDIR)
