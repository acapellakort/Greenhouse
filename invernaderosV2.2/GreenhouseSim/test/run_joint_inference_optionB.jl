# =============================================================================
# run_joint_inference_optionB.jl -- JOINT climate + growth t-walk (Note 25, Option B)
# =============================================================================
# Each MCMC step:
#   theta_v -> Vanthoor forward_map (days 0..N_DAYS-1, 5-min grid)
#           -> climate likelihood: hourly-binned T_air (+ RH) vs measured
#           -> growth model driven by VANTHOOR T_air at every 5-min sample
#              (CO2 and canopy light stay MEASURED -- Phase 1 coupling, T only)
#           -> growth likelihood: weekly yield, leaf number, derived LAI
#
# Same growth machinery as calibrate_growth.jl (simulate_growth_measured), same
# Vanthoor input construction as rmse_vanthoor_diagnostic_v3.jl. No src/ changes.
#
# Run from the GreenhouseSim folder:
#   julia --project=. test/run_joint_inference_optionB.jl smoke     # ~1-3 min, ALWAYS FIRST
#   SEED=1 CHAIN=150000 BURN=60000 julia --project=. test/run_joint_inference_optionB.jl full
#
# Parameters (13):
#   climate  nu4, gamma4, k_ground, U8_min_leak           (cover node ON from JSON)
#   growth   node_rate, veg_sink_max, LAI_max, set_rate   (J_max, Wf_max, ssd=15 fixed)
#   noise    sigma_T, sigma_RH, sigma_y, sigma_leaf, sigma_lai
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, JTwalk, JLD2, Plots, Printf, Statistics, Random

const MODE = get(ENV, "OPTB_MODE", isempty(ARGS) ? "smoke" : ARGS[1])   # "setup" = build inputs only (used by postprocess)
MODE in ("smoke", "full", "setup") || error("usage: ... run_joint_inference_optionB.jl [smoke|full]")
const SEED = parse(Int, get(ENV, "SEED", "1"))
Random.seed!(SEED)

# ── user knobs ───────────────────────────────────────────────────────────────
const PIPE2 = get(ENV, "PIPE2", "0") == "1"; const PIPE2_R = parse(Float64, get(ENV, "PIPE2_R", "0.6"))
const TEAM          = get(ENV, "TEAM", "AiCU")   # Note 54: other AGC-2018 compartments (Croperators, DeepGreens, iGrow, Sonoma, Reference(Growers))
const DAY0_SERIAL   = 43326.0
const N_DAYS        = parse(Int, get(ENV, "N_DAYS", "42"))   # 42 = AiCU window; 115 = full season (~3x cost)
const CHAIN_LENGTH  = parse(Int, get(ENV, "CHAIN", "200000"))
const BURN_IN       = parse(Int, get(ENV, "BURN", "80000"))
const SMOKE_STEPS   = 300
const CLIM_AGG_H    = 1.0         # climate likelihood on CLIM_AGG_H-hour bin means
const USE_RH        = true        # RH identifies gamma4 (as in the climate-only scripts)
const EXCLUDE_DAYS  = Int[]       # e.g. [24, 30, 41]: known structural spike days (Note 28/29)
const LAI_FIX       = 2.0         # Vanthoor LAI input (growth LAI NOT fed back -- Phase 1)

# Vanthoor input construction: identical to rmse_vanthoor_diagnostic_v3.jl best config
const PIPE_DELTA_K    = 25.0
const PIPE_MIN_VALID  = 5.0
const PIPE_WARM_TOUT  = 22.0
const CO2_AMB_PPM     = 390.0
const CO2_THR_PPM     = 50.0
const CO2_INJECT_FRAC = 0.3

# ── paths ────────────────────────────────────────────────────────────────────
const INV     = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DATA    = joinpath(INV, "invernaderos2.0", "data")
const GH_CSV  = joinpath(DATA, TEAM, "Greenhouse_climate.csv")
const MET_CSV = joinpath(DATA, "meteo.csv")
const VIP_CSV = joinpath(DATA, TEAM, "vip.csv")
const PROD_CSV = joinpath(DATA, TEAM, "Production.csv")
const CM_CSV  = joinpath(DATA, TEAM, "CropManagement.csv")
const CLM_JS  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
const CRP_JS  = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
const GRW_JS  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")
const N21_CHAIN = joinpath(@__DIR__, "posteriors", "AiCU_growth_chain.jld2")
const FIGDIR  = joinpath(@__DIR__, "figures"); mkpath(FIGDIR)
const POSTDIR = joinpath(@__DIR__, "posteriors"); mkpath(POSTDIR)

# ── parameters ───────────────────────────────────────────────────────────────
const QoI_V   = ["nu4", "gamma4", "k_ground", "U8_min_leak"]
const QoI_G   = ["node_rate", "veg_sink_max", "LAI_max", "set_rate"]
const QoI_S   = ["sigma_T", "sigma_RH", "sigma_y", "sigma_leaf", "sigma_lai"]
const QoI     = vcat(QoI_V, QoI_G, QoI_S)
const IV      = 1:4
const IG      = 5:8
const iσT, iσRH, iσy, iσl, iσL = 9, 10, 11, 12, 13

# [min, centre, max]. Climate priors = run_inference_venting_aicu.jl (+ U8_min_leak, Note 29);
# growth/noise priors = calibrate_growth.jl, so the joint posterior is comparable to Note 21.
const QoI_dict = Dict(
    "nu4"          => [1.0e-5, 4.46e-5, 1.0e-1],
    "gamma4"       => [20.0,   44.0,    300.0],
    "k_ground"     => [0.5,    24.5,    50.0],
    "U8_min_leak"  => [0.0,    0.05,    0.15],
    "node_rate"    => [0.02,   0.0919,  0.25],
    "veg_sink_max" => [2.0,    16.47,   40.0],
    "LAI_max"      => [1.5,    2.517,   4.5],
    "set_rate"     => [0.1,    1.302,   5.0],
    "sigma_T"      => [0.2,    2.5,     8.0],
    "sigma_RH"     => [0.5,    8.0,     30.0],
    "sigma_y"      => [0.02,   0.872,   3.0],
    "sigma_leaf"   => [1.0,    10.0,    40.0],   # lower bound as in calibrate_growth.jl
    "sigma_lai"    => [0.05,   0.5,     2.0],
)

# ── helpers (as in rmse_vanthoor_diagnostic_v3.jl) ───────────────────────────
const EXCEL_UNIX0 = 25569.0
const SIGMA_SB    = 5.670374419e-8
serial2unix(s) = (s - EXCEL_UNIX0) * 86400.0
_num(x) = (x isa Number && !(x isa Missing) && !isnan(Float64(x))) ? Float64(x) : NaN
stepfun(ts, vs) = t -> @inbounds vs[clamp(searchsortedlast(ts, t), 1, length(vs))]
_q2(T) = T > 0 ? 611.21*exp((18.678-T/234.5)*(T/(257.14+T))) :
                  611.21*exp((23.036-T/333.7)*(T/(279.82+T)))
_vpout(Tc, RH) = _q2(Tc) * clamp(RH, 1.0, 100.0) / 100.0

params = load_params(CLM_JS, CRP_JS)
@printf("Engine flags from JSON: use_cover_node=%.1f  h_cover_out=%.1f  (U8_min_leak in JSON=%.3f, overridden by theta_v)\n",
        params.use_cover_node, params.h_cover_out, params.U8_min_leak)
params.use_cover_node > 0.5 || @warn "use_cover_node is OFF in constants_climate.json -- Note 29 config expects 1"

# ── read data: Vanthoor inputs + climate obs + growth samples, one pass ─────
gh  = DataFrame(CSV.File(GH_CSV)); mt = DataFrame(CSV.File(MET_CSV)); vip = DataFrame(CSV.File(VIP_CSV))
# vip columns: AiCU has per-zone columns (X_1 .. X_6); the other teams have a single column X (Note 54)
function vcols(base, zs)
    nm = names(vip); z = [Symbol("$(base)_$(k)") for k in zs if "$(base)_$(k)" in nm]
    isempty(z) ? (base in nm ? [Symbol(base)] : Symbol[]) : z
end
const VC_HEAT = vcols("HeatTemp_Vip", (1, 2, 3, 5, 6))
n   = min(nrow(gh), nrow(mt), nrow(vip))
const win_lo = DAY0_SERIAL; const win_hi = DAY0_SERIAL + N_DAYS
const trans  = (1 - params.eta1) * params.tau1 * params.eta2

ts = Float64[]; U1v = Float64[]; U8v = Float64[]; U10v = Float64[]; U11v = Float64[]; U12v = Float64[]
Tpipev = Float64[]; Ig_v = Float64[]; Tout_v = Float64[]; Wind_v = Float64[]; VPout_v = Float64[]; Tsky_v = Float64[]
Tair_obs = Float64[]; CO2_obs = Float64[]; RH_obs = Float64[]
# growth samples (same validity rule as load_measured_climate: Tair & CO2 present)
gs_day = Int[]; gs_rel = Float64[]; gs_T = Float64[]; gs_C = Float64[]; gs_I = Float64[]

for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    To = _num(mt.Tout[i]); Ig = _num(mt.Iglob[i]); Ws = _num(mt.Windsp[i])
    Rho = _num(mt.Rhout[i]); Py = _num(mt.Pyrgeo[i])
    ToK = (isnan(To) ? 15.0 : To) + 273.15
    Igv = isnan(Ig) ? 0.0 : max(Ig, 0.0)
    Ldn = SIGMA_SB * ToK^4 + (isnan(Py) ? 0.0 : Py)
    Tskv = isnan(Py) ? (-0.4 + 273.15) : clamp((Ldn / SIGMA_SB)^0.25, ToK - 30.0, ToK + 2.0)

    vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i]); en = _num(gh.EnScr[i])
    al = _num(gh.AssimLight[i]); pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    U8  = clamp(((isnan(vl) ? 0.0 : vl) + (isnan(vw) ? 0.0 : vw)) / 200.0, 0.0, 1.0)
    U1  = isnan(en) ? 0.0 : clamp(en / 100.0, 0.0, 1.0)
    U12 = isnan(al) ? 0.0 : clamp(al / 100.0, 0.0, 1.0)

    hvals = filter(!isnan, [_num(vip[i, c]) for c in VC_HEAT])
    T_heat_C = isempty(hvals) ? 20.0 : mean(hvals)

    pg_c = (!isnan(pg) && pg > PIPE_MIN_VALID) ? pg : -Inf
    pl_c = (!isnan(pl) && pl > PIPE_MIN_VALID) ? pl : -Inf
    pipe_max = max(pg_c, pl_c)
    # Note 54 (default-off, ENV PIPE2=1): when both circuits are warm, replace max(PipeLow, PipeGrow) by the single
    # pipe temperature that delivers the same convective heat as both, h ~ dT^1.32 (dT vs the LOGGED air),
    # the grow pipe weighted by PIPE2_R (heat per K^1.32 relative to the rail pipe; default 0.6).
    if PIPE2 && pg_c > -Inf && pl_c > -Inf
        Tal = _num(gh.Tair[i])
        if !isnan(Tal) && pg_c > Tal && pl_c > Tal
            pipe_max = Tal + ((pl_c - Tal)^1.32 + PIPE2_R * (pg_c - Tal)^1.32)^(1 / 1.32)
        end
    end
    To_val = isnan(To) ? 15.0 : To
    Tpipe = pipe_max > -Inf ? pipe_max + 273.15 :
            (To_val < T_heat_C - 2.0 && To_val < PIPE_WARM_TOUT) ? T_heat_C + PIPE_DELTA_K + 273.15 :
            To_val + 2.0 + 273.15

    co2vip = _num(vip.CO2_Vip[i])
    U10 = (isnan(co2vip) || co2vip <= CO2_AMB_PPM + CO2_THR_PPM) ? 0.0 : CO2_INJECT_FRAC

    Ta = _num(gh.Tair[i]); Co = _num(gh.CO2air[i]); Rh = _num(gh.RHair[i])
    push!(ts, serial2unix(g))
    push!(Ig_v, Igv); push!(Tout_v, ToK); push!(Wind_v, isnan(Ws) ? 2.0 : Ws)
    push!(VPout_v, _vpout(To_val, isnan(Rho) ? 80.0 : Rho)); push!(Tsky_v, Tskv)
    push!(U1v, U1); push!(U8v, U8); push!(U10v, U10); push!(U11v, T_heat_C + 273.15)
    push!(U12v, U12); push!(Tpipev, Tpipe)
    push!(Tair_obs, isnan(Ta) ? NaN : Ta + 273.15)
    push!(CO2_obs, isnan(Co) ? NaN : Co * 536.4 / ((isnan(Ta) ? 20.0 : Ta) + 273.15))
    push!(RH_obs, isnan(Rh) ? NaN : Rh)

    if !isnan(Ta) && !isnan(Co)
        push!(gs_day, Int(floor(g - DAY0_SERIAL)))
        push!(gs_rel, serial2unix(g) - serial2unix(DAY0_SERIAL))
        push!(gs_T, Ta + 273.15)
        push!(gs_C, Co * 536.4 / (Ta + 273.15))
        push!(gs_I, trans * Igv + params.alpha12 * ((isnan(al) ? 0.0 : al) / 100.0))
    end
end
median(Wind_v) > 25.0 && (Wind_v ./= 3.6)

const t0 = ts[1]
abs(t0 - serial2unix(DAY0_SERIAL)) < 1.0 || error("first window row is not at DAY0 midnight; grid mapping assumes it is")
trel = ts .- t0
A(v) = stepfun(ts, v); Z(v) = stepfun(trel, v); zc = Z(zeros(length(trel)))

weather  = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, _t -> 18.0 + 273.15,
                         A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), _t -> LAI_FIX)
controls = ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8v), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipev))
const base = SimBase(params, weather, controls, t0)

i0  = findfirst(!isnan, Tair_obs)
Ta0 = Tair_obs[i0]; Co0 = isnan(CO2_obs[i0]) ? 575.0 : CO2_obs[i0]; Rh0 = isnan(RH_obs[i0]) ? 80.0 : RH_obs[i0]
const u0 = [Ta0, Ta0, Co0, (Rh0 / 100.0) * Pws(Ta0), 0.0, 18.0 + 273.15]

# ── prediction grid (same layout as forward_map/_run) ────────────────────────
const days  = collect(0.0:(N_DAYS - 1))
const tspan = (0.0, 86400.0)
const teval = 0.0:300.0:86400.0
const NPD   = length(teval)                       # 289 saved points per day
const NCOL  = 1 + N_DAYS * NPD                    # +1 = initial state column
col_of(d, j) = 2 + d * NPD + j                    # day d (0-based), grid point j (0..288)

predtimes = vcat(t0, [t0 + d * 86400.0 + tj for d in days for tj in teval])
Tobs_g  = A(Tair_obs).(predtimes)
RHobs_g = A(RH_obs).(predtimes)

# hour bins for the climate likelihood (bins need >= half their points observed)
function make_bins(obs)
    nb = round(Int, 24 / CLIM_AGG_H)
    bins = Vector{Vector{Int}}(); means = Float64[]
    for d in 0:(N_DAYS - 1)
        d in EXCLUDE_DAYS && continue
        for b in 0:(nb - 1)
            cols = [col_of(d, j) for j in 0:(NPD - 1) if min(floor(Int, j * 300 / (3600 * CLIM_AGG_H)), nb - 1) == b]
            ok = [c for c in cols if !isnan(obs[c])]
            length(ok) >= length(cols) / 2 || continue
            push!(bins, ok); push!(means, mean(obs[ok]))
        end
    end
    return bins, means
end
const binsT, TbarObs   = make_bins(Tobs_g)
const binsRH, RHbarObs = make_bins(RHobs_g)
binmean(v, ix) = (s = 0.0; @inbounds for i in ix; s += v[i]; end; s / length(ix))

# ── growth season driven by measured CO2/light, T overwritten per step ───────
samp_meas = Vector{Matrix{Float64}}(undef, N_DAYS)
samp_col  = Vector{Vector{Int}}(undef, N_DAYS)
for d in 0:(N_DAYS - 1)
    ix = findall(==(d), gs_day)
    M = Matrix{Float64}(undef, 3, length(ix)); c = Vector{Int}(undef, length(ix))
    for (k, i) in enumerate(ix)
        M[1, k] = gs_T[i]; M[2, k] = gs_C[i]; M[3, k] = gs_I[i]
        j = clamp(round(Int, (gs_rel[i] - 86400.0 * d) / 300.0), 0, NPD - 1)
        c[k] = col_of(d, j)
    end
    samp_meas[d + 1] = M; samp_col[d + 1] = c
end
const work      = [copy(M) for M in samp_meas]          # row 1 overwritten each step
const season_B  = MeasuredSeason(N_DAYS, work, 300.0)    # holds `work` by reference
const season_A  = MeasuredSeason(N_DAYS, samp_meas, 300.0)

gp_base = load_growth_params(GRW_JS)
gp_base = update_params(gp_base, ["set_start_day"], [15.0])

# growth observations restricted to the window
days_obs, y_obs = load_production(PROD_CSV)
cum_at(dd, cum, d) = (d <= dd[1] ? float(cum[1]) : d >= dd[end] ? float(cum[end]) :
    let i = searchsortedlast(dd, d)
        dd[i] == d ? float(cum[i]) : cum[i] + (d - dd[i]) / (dd[i+1] - dd[i]) * (cum[i+1] - cum[i])
    end)
weekly_incr(dd, cum, e) = diff([cum_at(dd, cum, x) for x in e])
const LAST_DAY = min(N_DAYS - 1, days_obs[end])
const edges = let e = collect(0:7:LAST_DAY); e[end] < LAST_DAY && push!(e, LAST_DAY); e end
const Δy_obs = weekly_incr(days_obs, y_obs, edges)
cm_days_all, cm_N_all, cm_L_all = load_cropmanagement(CM_CSV)
keep = cm_days_all .<= N_DAYS - 1
const cm_days = cm_days_all[keep]; const cm_N = cm_N_all[keep]; const cm_L = cm_L_all[keep]
const Nwk = length(Δy_obs); const Ncm = length(cm_days)
at_day(v, d) = v[min(d, length(v) - 1) + 1]

# ── energy ───────────────────────────────────────────────────────────────────
PriorSupp(z) = all(QoI_dict[QoI[k]][1] < z[k] < QoI_dict[QoI[k]][3] for k in eachindex(z))

function growth_terms(res, z)
    σy, σl, σL = z[iσy], z[iσl], z[iσL]
    ry = Δy_obs .- weekly_incr(res.day, res.yield_FW, edges)
    rl = cm_N .- [at_day(res.node, d) for d in cm_days]
    rL = cm_L .- [at_day(res.LAI,  d) for d in cm_days]
    Uy = Nwk * log(σy) + sum(abs2, ry) / (2σy^2)
    Ul = Ncm * log(σl) + sum(abs2, rl) / (2σl^2)
    UL = Ncm * log(σL) + sum(abs2, rL) / (2σL^2)
    return Uy, Ul, UL
end

const BAD = (U = Inf, UT = Inf, URH = Inf, Uy = Inf, Ul = Inf, UL = Inf)

function energy_parts(z)
    local T2, RH
    try
        _, T2, RH, _ = forward_map(z[IV], QoI_V, base, days, tspan, teval; u0 = copy(u0))
    catch
        return BAD
    end
    (length(T2) == NCOL && all(isfinite, T2) && all(isfinite, RH)) || return BAD

    σT = z[iσT]; sT = 0.0
    for k in eachindex(binsT); sT += (binmean(T2, binsT[k]) - TbarObs[k])^2; end
    UT = length(binsT) * log(σT) + sT / (2σT^2)
    URH = 0.0
    if USE_RH
        σR = z[iσRH]; sR = 0.0
        for k in eachindex(binsRH); sR += (binmean(RH, binsRH[k]) - RHbarObs[k])^2; end
        URH = length(binsRH) * log(σR) + sR / (2σR^2)
    end

    # Option B coupling: Vanthoor T_air replaces measured T at every growth sample
    for d in 1:N_DAYS
        M = work[d]; c = samp_col[d]
        @inbounds for k in eachindex(c); M[1, k] = T2[c[k]]; end
    end
    gp  = update_params(gp_base, QoI_G, z[IG])
    res = try simulate_growth_measured(params, gp, season_B; LAI0 = 0.5) catch; return BAD end
    Uy, Ul, UL = growth_terms(res, z)
    U = UT + URH + Uy + Ul + UL
    return isfinite(U) ? (; U, UT, URH, Uy, Ul, UL) : BAD
end
energy(z) = energy_parts(z).U

# ── progress bar + checkpoints (JTwalk untouched) ────────────────────────────
# JTwalk.Run! calls Supp exactly once per iteration (every kernel with the default
# weights), plus twice in SetupInitialValues! for x0/xp0. So a counting Supp wrapper
# gives the iteration number. The same hook periodically saves the chain rows written
# so far (obj.Output[1:it, :]) so a crash does not lose a long run.
mutable struct Progress
    total::Int; calls::Int; t0::Float64; tlast::Float64; tckpt::Float64
    Umin::Float64; ckpt::String
end
const PROG   = Progress(0, 0, 0.0, 0.0, 0.0, Inf, "")
const OBJREF = Ref{Any}(nothing)
const IS_TTY = stdout isa Base.TTY
const PROGRESS_SEC = parse(Float64, get(ENV, "PROGRESS_SEC", MODE == "smoke" ? "10" : (IS_TTY ? "30" : "300")))
const CKPT_MIN     = parse(Float64, get(ENV, "CKPT_MIN", "30"))     # 0 disables checkpoints

fmt_t(s) = !isfinite(s) ? "?" : s < 3600 ? @sprintf("%dm%02ds", floor(Int, s / 60), floor(Int, mod(s, 60))) :
                                          @sprintf("%dh%02dm", floor(Int, s / 3600), floor(Int, mod(s, 3600) / 60))
current_iter() = clamp(PROG.calls - 2, 0, PROG.total)

function start_progress!(obj, total, ckpt)
    OBJREF[] = obj; PROG.total = total; PROG.calls = 0; PROG.Umin = Inf; PROG.ckpt = ckpt
    PROG.t0 = PROG.tlast = PROG.tckpt = time()
end

function report!(force = false)
    now = time()
    (force || now - PROG.tlast >= PROGRESS_SEC) || return
    PROG.tlast = now
    it = current_iter(); el = now - PROG.t0
    rate = it > 0 ? el / it : NaN; eta = it > 0 ? rate * (PROG.total - it) : NaN
    frac = it / PROG.total; nb = round(Int, 30 * frac)
    line = @sprintf("[%s%s] %5.1f%%  %d/%d  elapsed %s  %.0f ms/it  ETA %s  Umin %.1f",
                    "#"^nb, "-"^(30 - nb), 100 * frac, it, PROG.total, fmt_t(el), 1e3 * rate, fmt_t(eta), PROG.Umin)
    IS_TTY ? print("\r", line) : println(line)
    flush(stdout)
end

function checkpoint!()
    obj = OBJREF[]; (obj === nothing || isempty(PROG.ckpt)) && return
    rows = min(current_iter(), size(obj.Output, 1))      # rows 1..it are filled
    rows < 2 && return
    tmp = PROG.ckpt * ".tmp"
    JLD2.jldopen(tmp, "w") do f
        f["chain"] = obj.Output[1:rows, :]; f["QoI"] = QoI; f["iter"] = rows - 1
        f["burn_in"] = BURN_IN; f["chain_length"] = CHAIN_LENGTH; f["n_days"] = N_DAYS; f["seed"] = SEED
    end
    mv(tmp, PROG.ckpt; force = true)
end

function tracked_supp(z)
    PROG.calls += 1
    report!()
    if CKPT_MIN > 0 && time() - PROG.tckpt >= 60 * CKPT_MIN
        PROG.tckpt = time(); checkpoint!()
    end
    return PriorSupp(z)
end
tracked_energy(z) = (U = energy(z); U < PROG.Umin && (PROG.Umin = U); U)

# ── warm start: Note 21 posterior mean (growth + growth sigmas), Note 28/29 (climate) ──
centre = [QoI_dict[q][2] for q in QoI]
if isfile(N21_CHAIN)
    c21 = JLD2.load(N21_CHAIN)
    m21 = vec(mean(c21["chain"][c21["burn_in"]:end, 1:end-1], dims = 1))
    for (k, q) in enumerate(c21["QoI"])
        i = findfirst(==(q), QoI); isnothing(i) || (centre[i] = m21[k])
    end
    println("Warm start from Note 21 chain: ", N21_CHAIN)
else
    @warn "Note 21 chain not found; warm start uses QoI_dict centres"
end
function jitter(c)
    for _ in 1:1000
        x = c .* (1 .+ 0.02 .* randn(length(c)))
        PriorSupp(x) && return x
    end
    error("could not draw a start point inside the prior")
end

println("\n", "="^70)
@printf("Option B joint inference -- %s, days 0-%d, mode=%s, seed=%d\n", TEAM, N_DAYS - 1, MODE, SEED)
println("="^70)
@printf("Climate bins: T %d, RH %d  (%.0f-h means; excluded days %s)\n",
        length(binsT), length(binsRH), CLIM_AGG_H, string(EXCLUDE_DAYS))
@printf("Growth obs:   %d weekly yield bins (%d non-zero), %d leaf + %d LAI points\n",
        Nwk, count(>(1e-9), Δy_obs), Ncm, Ncm)
@printf("Growth samples per day: min %d, median %d\n",
        minimum(size.(samp_meas, 2)), round(Int, median(size.(samp_meas, 2))))

# ═════════════════════════════════════════════════════════════════════════════
if MODE == "smoke"
    # 1. our hand-built season must reproduce load_measured_climate exactly
    lib = load_measured_climate(GH_CSV, MET_CSV, params)
    lib_cut = MeasuredSeason(N_DAYS, lib.day_samples[1:N_DAYS], lib.dt)
    same_n  = all(size(lib.day_samples[d], 2) == size(samp_meas[d], 2) for d in 1:N_DAYS)
    maxdiff = same_n ? maximum(maximum(abs.(lib.day_samples[d] .- samp_meas[d]); init = 0.0) for d in 1:N_DAYS) : NaN
    gp_c = update_params(gp_base, QoI_G, centre[IG])
    rA = simulate_growth_measured(params, gp_c, season_A; LAI0 = 0.5)
    rL = simulate_growth_measured(params, gp_c, lib_cut;  LAI0 = 0.5)
    @printf("\n[check 1] season vs load_measured_climate: same sample counts=%s  max|diff|=%.2e  yield diff=%.2e\n",
            same_n, maxdiff, abs(rA.yield_FW[end] - rL.yield_FW[end]))

    # 2. energy at the warm-start centre (first call compiles; time the second)
    energy_parts(centre)
    t_e = @elapsed (ep = energy_parts(centre))
    @printf("\n[check 2] energy at warm start: U=%.2f  (UT=%.2f URH=%.2f Uy=%.2f Uleaf=%.2f ULAI=%.2f)\n",
            ep.U, ep.UT, ep.URH, ep.Uy, ep.Ul, ep.UL)
    @printf("          one energy evaluation: %.1f ms\n", 1e3 * t_e)

    # 3. Option A vs Option B at the same growth parameters (the size of the coupling effect)
    rB = simulate_growth_measured(params, gp_c, season_B; LAI0 = 0.5)   # season_B now holds Vanthoor T
    UyA, UlA, ULA = growth_terms(rA, centre)
    @printf("\n[check 3] growth at Note 21 posterior mean, day %d:\n", N_DAYS - 1)
    @printf("          measured T (A): yield %.2f kg/m2, nodes %.1f, U_growth %.2f\n",
            rA.yield_FW[end], rA.node[end], UyA + UlA + ULA)
    @printf("          Vanthoor T (B): yield %.2f kg/m2, nodes %.1f, U_growth %.2f\n",
            rB.yield_FW[end], rB.node[end], ep.Uy + ep.Ul + ep.UL)
    @printf("          observed:       yield %.2f kg/m2 (day %d)\n", cum_at(days_obs, y_obs, LAST_DAY), LAST_DAY)
    mTv = mean(mean(work[d][1, :]) for d in 1:N_DAYS if size(work[d], 2) > 0)
    mTm = mean(mean(samp_meas[d][1, :]) for d in 1:N_DAYS if size(samp_meas[d], 2) > 0)
    @printf("          mean T at growth samples: Vanthoor %.2f C vs measured %.2f C (bias %+.2f K)\n",
            mTv - 273.15, mTm - 273.15, mTv - mTm)

    # 4. short t-walk: acceptance and projected wall time
    obj = jtwalk(n = length(QoI), U = tracked_energy, Supp = tracked_supp)
    x0s, xp0s = jitter(centre), jitter(centre)
    start_progress!(obj, SMOKE_STEPS, "")                 # no checkpoint in smoke
    t_w = @elapsed Run!(obj, T = SMOKE_STEPS, x0 = x0s, xp0 = xp0s)
    report!(true); println()
    PROG.calls - 2 == SMOKE_STEPS || @warn "Supp-call count ($(PROG.calls - 2)) != steps ($SMOKE_STEPS): progress counter is off"
    O = obj.Output
    acc = mean(any(O[k+1, 1:end-1] .!= O[k, 1:end-1]) for k in 1:size(O, 1)-1)
    @printf("\n[check 4] %d t-walk steps in %.1f s  (%.1f ms/step), moved-fraction %.2f\n",
            SMOKE_STEPS, t_w, 1e3 * t_w / SMOKE_STEPS, acc)
    @printf("          energy start %.1f -> end %.1f\n", O[1, end], O[end, end])
    @printf("          projected full run (%d steps): %.1f h\n", CHAIN_LENGTH, t_w / SMOKE_STEPS * CHAIN_LENGTH / 3600)
    println("\nSmoke done. Paste this output before launching `full`.")
    exit(0)
end

# ═════════════════════════════════════════════════════════════════════════════
# FULL RUN
if MODE == "full"
tag = "$(TEAM)_joint_optionB_d$(N_DAYS)_seed$(SEED)"
obj = jtwalk(n = length(QoI), U = tracked_energy, Supp = tracked_supp)
println("Inferring: ", QoI, "  chain=", CHAIN_LENGTH)
@printf("Progress every %.0f s; checkpoint every %.0f min -> %s\n", PROGRESS_SEC, CKPT_MIN,
        joinpath(POSTDIR, tag * "_ckpt.jld2"))
x0f, xp0f = jitter(centre), jitter(centre)
start_progress!(obj, CHAIN_LENGTH, joinpath(POSTDIR, tag * "_ckpt.jld2"))
@time Run!(obj, T = CHAIN_LENGTH, x0 = x0f, xp0 = xp0f)
report!(true); println()
JLD2.@save joinpath(POSTDIR, tag * ".jld2") chain=obj.Output QoI=QoI burn_in=BURN_IN chain_length=CHAIN_LENGTH n_days=N_DAYS seed=SEED clim_agg_h=CLIM_AGG_H use_rh=USE_RH exclude_days=EXCLUDE_DAYS
println("Chain saved: ", joinpath(POSTDIR, tag * ".jld2"))

post = obj.Output[BURN_IN:end, 1:end-1]
pm = vec(mean(post, dims = 1)); ps = vec(std(post, dims = 1))
sd21 = Dict{String,Float64}(); mean21 = Dict{String,Float64}()
if isfile(N21_CHAIN)
    c21 = JLD2.load(N21_CHAIN); p21 = c21["chain"][c21["burn_in"]:end, 1:end-1]
    for (k, q) in enumerate(c21["QoI"]); mean21[q] = mean(p21[:, k]); sd21[q] = std(p21[:, k]); end
end
println("\n--- joint posterior (mean ± sd) ---   [shift vs Note 21 in Note-21 SDs]")
for (i, q) in enumerate(QoI)
    s = haskey(sd21, q) ? @sprintf("   [%+.1f sd]", (pm[i] - mean21[q]) / sd21[q]) : ""
    @printf("  %-13s = %.4g ± %.3g%s\n", q, pm[i], ps[i], s)
end
acc = mean(any(obj.Output[k+1, 1:end-1] .!= obj.Output[k, 1:end-1]) for k in 1:size(obj.Output, 1)-1)
@printf("moved-fraction: %.3f\n", acc)

ep = energy_parts(pm)
_, T2m, _, _ = forward_map(pm[IV], QoI_V, base, days, tspan, teval; u0 = copy(u0))
mk = .!isnan.(Tobs_g); mk[1] = false
@printf("At posterior mean: RMSE_T (5-min) %.2f K, bias %+.2f K; U parts: T %.1f RH %.1f y %.1f leaf %.1f LAI %.1f\n",
        sqrt(mean((T2m[mk] .- Tobs_g[mk]).^2)), mean(T2m[mk] .- Tobs_g[mk]), ep.UT, ep.URH, ep.Uy, ep.Ul, ep.UL)
res = simulate_growth_measured(params, update_params(gp_base, QoI_G, pm[IG]), season_B; LAI0 = 0.5)
@printf("Yield day %d: model %.2f vs obs %.2f kg/m2\n", LAST_DAY, res.yield_FW[LAST_DAY + 1], cum_at(days_obs, y_obs, LAST_DAY))

pE = plot(obj.Output[:, end], label = "", title = "Joint energy", xlabel = "iter"); vline!(pE, [BURN_IN], ls = :dash, lc = :red, label = "")
savefig(pE, joinpath(FIGDIR, tag * "_energy.png"))
for (i, q) in enumerate(QoI)
    p1 = plot(obj.Output[:, i], label = "", title = "chain: $q"); vline!(p1, [BURN_IN], ls = :dash, lc = :red, label = "")
    savefig(p1, joinpath(FIGDIR, tag * "_$(q)_chain.png"))
end
hrs = (predtimes .- t0) ./ 86400
pT = plot(hrs, Tobs_g .- 273.15, lc = :black, lw = 2, label = "measured", ylabel = "T air (C)", xlabel = "day", size = (900, 350))
plot!(pT, hrs, collect(T2m) .- 273.15, lc = :darkorange, label = "joint posterior mean")
savefig(pT, joinpath(FIGDIR, tag * "_Tfit.png"))
pW = scatter(edges[2:end], Δy_obs, label = "observed", xlabel = "day", ylabel = "kg FW/m2/wk")
plot!(pW, edges[2:end], weekly_incr(res.day, res.yield_FW, edges), lw = 2, label = "model (Vanthoor T)")
savefig(pW, joinpath(FIGDIR, tag * "_yield.png"))
println("Figures: ", FIGDIR, "/", tag, "_*.png")
end  # MODE == "full"
