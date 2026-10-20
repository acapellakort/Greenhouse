# =============================================================================
# run_climate_calibration_cl.jl -- CLOSED-LOOP climate calibration (Note 35 -> Note 36)
# =============================================================================
# Run from the GreenhouseSim folder:
#   DESIGN=blocks julia --project=. test/run_climate_calibration_cl.jl smoke          # ALWAYS FIRST
#   DESIGN=blocks SEED=1 CHAIN=100000 BURN=40000 RELTOL=1e-3 \
#       julia --project=. test/run_climate_calibration_cl.jl full                       # only after OK
#
# Same structure as run_climate_calibration.jl (Note 33), with three changes:
#   1. vent_ctrl = 1: the roof vents are a model STATE (u[7]) driven by the grower's vent
#      setpoint and humidity-deficit setpoint (ControlInputsX), no longer replayed from the log.
#   2. The LOGGED vent opening becomes a second observable: hourly bin means, own noise sigma_U.
#   3. Each day is solved here with an explicit `reltol` (env RELTOL, default 1e-3), because
#      GreenhouseSim._run takes no tolerance argument. The smoke mode checks 1e-3 vs 1e-4 vs 1e-5.
#
# Engine flags fixed ON : heat_ctrl (Kp 10, Ki 10/1800), vent_screen, stomata_fix, vapour_fix,
#                         cond_heat with f_cond = 0 (i.e. latent heat off, as Note 33 found ~0),
#                         vent_ctrl. U8_min_leak = 0, nu4 = 4.46e-5, LAI(t) input, measured pipe.
# Engine flags fixed OFF: pipe_lag, black_scr (Note 35: no gain at tested values).
# Free (12): gamma4, vent_exp, vent_scr_perm, hc_dT, k_ground,
#            vc_gain, vc_e0, vc_tau, vc_hum,  sigma_T, sigma_RH, sigma_U.
# DESIGN: blocks only by default (TRAIN 7-13, 28-34, 63-69, 91-97; VAL 42-48, 77-83; TEST rest);
#         early is accepted for comparison. COMMON TEST = 49-62, 70-76, 84-90, 98-114.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
get(ENV, "PRIORS", "") in ("n54", "n58", "n61", "n61b", "n79") && !haskey(ENV, "PIPE2") && (ENV["PIPE2"] = "1")   # Note 54: both heating circuits
ENV["N_DAYS"]    = "115"
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))   # season inputs (ts, U*v, Tair_obs, ...)

const CMODE  = get(ENV, "CL_MODE", isempty(ARGS) ? "smoke" : ARGS[1])   # "setup" = build inputs only (used by diagnostics)
CMODE in ("smoke", "full", "setup") || error("usage: ... run_climate_calibration_cl.jl [smoke|full]")
const DESIGN = get(ENV, "DESIGN", "blocks")
DESIGN in ("blocks", "early") || error("DESIGN must be blocks or early")
const CSEED  = parse(Int, get(ENV, "SEED", "1")); Random.seed!(CSEED)
const CCHAIN = parse(Int, get(ENV, "CHAIN", "100000"))
const CBURN  = parse(Int, get(ENV, "BURN", "40000"))
const RELTOL = parse(Float64, get(ENV, "RELTOL", "1e-3"))
const CSMOKE = 200
const PRIORS = get(ENV, "PRIORS", "n36")     # n36 = Note 36 priors; n37 = Note 34/37-informed controller priors + free vc_hband
PRIORS in ("n36", "n37", "n39", "n43", "n46", "n49a", "n49b", "n54", "n58", "n61", "n61b", "n79") || error("PRIORS must be n36, n37, n39, n43, n46, n49a, n49b, n54, n58, n61, n61b or n79")
const N79 = PRIORS == "n79"   # Note 79: a-priori geometry fixed, cover condensation + crop respiration on, LAI and respiration from the n78 crop posterior, U10 from logged daily dosage, no bias terms   # n61b = Note 62: n61 with s_rad = 0 and f_mass up to 0.5 (physics vs sensor)
const N61P = PRIORS in ("n61", "n61b", "n79")   # n61 = Note 61: n58 + air_cap_mult + s_rad (measurement), TRAIN + block 110-114   # n58 = Note 58d: OPEN LOOP, logged lee/wind vents (vent_split), lamp_full + screen loss   # n54 = Note 54b: n49b + full lamp heat (lamp_full) + both pipe circuits (PIPE2)   # n49a/b = Note 50: physical caps, free vent wind term; a = with daytime observation offset dday   # n46 = Notes 46-48: SF structure + solar-absorbing internal mass   # n43 = Note 43 structure: logged pipe + PI gaps, deep-soil state, logged min vent   # n39 = n37 + free vc_min_night (Note 40)
const CTAG   = "$(TEAM)_climcal_cl_$(DESIGN)$(PRIORS == "n36" ? "" : "_" * PRIORS)_seed$(CSEED)"

for k in (:vent_ctrl, :vc_gain, :vc_e0, :vc_tau, :vc_hum, :vc_hband, :pipe_lag, :black_scr, :heat_ctrl, :cond_heat)
    hasproperty(params, k) || error("constants_climate.json has no $k entry (Note 34 flags missing?)")
end

wk(a, b) = collect(a:b)
const TRAIN_BLOCKS = DESIGN == "blocks" ? vcat([(7, 13), (28, 34), (63, 69), (91, 97)], N61P ? [(110, 114)] : Tuple{Int,Int}[]) :
                                          [(7, 13), (21, 27), (35, 41)]   # n61: season-end block added (Note 61, mode separation)
const VAL_BLOCKS   = DESIGN == "blocks" ? [(42, 48), (77, 83)] : [(28, 34)]
const TRAIN_DAYS   = reduce(vcat, [wk(a, b) for (a, b) in TRAIN_BLOCKS])
const VAL_DAYS     = reduce(vcat, [wk(a, b) for (a, b) in VAL_BLOCKS])
const TEST_DAYS    = DESIGN == "blocks" ? setdiff(0:114, TRAIN_DAYS, VAL_DAYS) : collect(42:114)
const COMMON_TEST  = N61P ? vcat(wk(49, 62), wk(70, 76), wk(84, 90), wk(98, 109)) :   # n61: 110-114 is TRAIN
                                         vcat(wk(49, 62), wk(70, 76), wk(84, 90), wk(98, 114))

# ── inputs (as diagnose_closed_loop.jl) ─────────────────────────────────────
# Note 86: AiCU's vip.csv holds the setpoints of ALL six compartments (_1 iGrow, _2 DeepGreens, _3 Reference, _4 AiCU,
# _5 Sonoma, _6 Croperators); earlier priors averaged them. n79 uses AiCU's own column _4 (its heating setpoint _4 is empty).
const VIPZ = N79 ? (4,) : 1:6
const VC_LEEW = vcols("VentLeew_Vip", VIPZ); const VC_WIND = vcols("VentWind_Vip", VIPZ); const VC_VMIN = vcols("VentLeewMin_Vip", VIPZ)
LeeZ = Float64[]; WindZ = Float64[]; gdayv = Int[]      # Note 58d: logged lee / wind openings (fractions, NaN -> 0)
VentMinv = Float64[]; srcv = Int[]; VentSPv = Float64[]; HumSPv = Float64[]; Blackv = Float64[]; U8log = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    To = _num(mt.Tout[i]); To_val = isnan(To) ? 15.0 : To
    pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    hv = filter(!isnan, [_num(vip[i, c]) for c in VC_HEAT])
    Th = isempty(hv) ? 20.0 : mean(hv)
    valid = (!isnan(pg) && pg > PIPE_MIN_VALID) || (!isnan(pl) && pl > PIPE_MIN_VALID)
    push!(gdayv, floor(Int, g))
    push!(srcv, valid ? 1 : (To_val < Th - 2.0 && To_val < PIPE_WARM_TOUT) ? 2 : 3)
    vs = filter(!isnan, vcat([_num(vip[i, c]) for c in VC_LEEW], [_num(vip[i, c]) for c in VC_WIND]))
    push!(VentSPv, (isempty(vs) ? 25.0 : mean(vs)) + 273.15)
    hd = _num(vip.HumDef_Vip[i]); push!(HumSPv, isnan(hd) ? 0.0 : hd)
    bs = _num(gh.BlackScr[i]); push!(Blackv, isnan(bs) ? 0.0 : clamp(bs / 100, 0.0, 1.0))
    vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i])
    push!(U8log, (isnan(vl) || isnan(vw)) ? NaN : clamp((vl + vw) / 200, 0.0, 1.0))
    push!(LeeZ, isnan(vl) ? 0.0 : clamp(vl / 100, 0.0, 1.0)); push!(WindZ, isnan(vw) ? 0.0 : clamp(vw / 100, 0.0, 1.0))
    vmn = filter(!isnan, [_num(vip[i, c]) for c in VC_VMIN])
    push!(VentMinv, isempty(vmn) ? 0.0 : clamp(mean(vmn) / 200, 0.0, 1.0))
end
length(srcv) == length(ts) || error("row filter mismatch")
Tpipe_meas = [srcv[i] == 2 ? Tout_v[i] + 2.0 : Tpipev[i] for i in eachindex(Tpipev)]
lai_d = vcat(0.0, Float64.(cm_days_all)); lai_v = vcat(0.5, cm_L_all)
lai_of_day(x) = x <= lai_d[1] ? lai_v[1] : x >= lai_d[end] ? lai_v[end] :
    (i = searchsortedlast(lai_d, x); lai_v[i] + (x - lai_d[i]) / (lai_d[i+1] - lai_d[i]) * (lai_v[i+1] - lai_v[i]))
LAI_row = [lai_of_day((t - t0) / 86400.0) for t in ts]
# Note 79: LAI and crop respiration from the n78 crop posterior (test/crop_posterior_series.jl); CO2 valve from the
# logged daily dosage: on-times from the CO2 setpoint rule, scaled per day so valve hours = CO2_dosage / 0.0128 kg m-2 h-1
CropRespv = zeros(length(ts)); U10_scale = Float64[]
if N79
    csf = joinpath(@__DIR__, "posteriors", "AiCU_crop_series_n78.csv")
    isfile(csf) || error("missing $csf -- run test/crop_posterior_series.jl first")
    cs = DataFrame(CSV.File(csf))
    cday(x) = clamp(floor(Int, x), cs.day[1], cs.day[end]) - cs.day[1] + 1
    LAI_row   = [cs.LAI_floor[cday((t - t0) / 86400.0)] for t in ts]
    CropRespv = [cs.CropResp20[cday((t - t0) / 86400.0)] for t in ts]
    rc = DataFrame(CSV.File(joinpath(DATA, TEAM, "ResourceCalculations.csv")))
    dos = Dict(floor(Int, rc.time[i]) => rc.CO2_dosage[i] for i in 1:nrow(rc) if !ismissing(rc.CO2_dosage[i]) && !isnan(rc.CO2_dosage[i]))
    dtr = length(ts) > 1 ? median(diff(ts)) / 3600 : 300 / 3600          # h per row
    for gd in unique(gdayv)
        rows = findall(==(gd), gdayv); hlog = get(dos, gd, NaN) / 0.0128
        on = [i for i in rows if U10v[i] > 0]; hon = length(on) * dtr
        if isnan(hlog)
            sc = 1.0; for i in on; U10v[i] = 1.0; end                          # no dosage logged: valve fully open when on
        elseif hon > 0
            sc = hlog / hon; for i in on; U10v[i] = min(sc, 1.0); end
        else
            sc = NaN; dayr = [i for i in rows if Ig_v[i] > 50.0]                # dosage but no on-signal: spread over daylight
            for i in dayr; U10v[i] = isempty(dayr) ? 0.0 : min(hlog / (length(dayr) * dtr), 1.0); end
        end
        push!(U10_scale, sc)
    end
    sfin = filter(isfinite, U10_scale)
    @printf("Note 79 inputs: LAI_floor %.2f-%.2f, CropResp20 mean %.4f mg/m2/s; CO2 valve scale per day median %.2f (%d of %d days capped at 1)\n",
            minimum(LAI_row), maximum(LAI_row), mean(CropRespv), isempty(sfin) ? NaN : median(sfin), count(>(1.0), sfin), length(U10_scale))
end
const weather_c = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, _t -> 18.0 + 273.15,
                                A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), A(LAI_row))
U8z = [isnan(v) ? 0.0 : v for v in U8log]
const N46P = PRIORS in ("n46", "n49a", "n49b", "n54", "n58", "n61", "n61b", "n79")   # Notes 46-48 structure
const OLP  = PRIORS in ("n58", "n61", "n61b", "n79")     # Note 58d: open-loop calibration (vents are logged inputs, no vent likelihood)
const N43P = PRIORS in ("n43", "n46", "n49a", "n49b", "n54", "n58", "n61", "n61b", "n79")     # structures with the Note 43 flags, VentMin input and deep-soil state
const ctrl_x = N79 ?
    ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, Z(WindZ), Z(LeeZ), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                   Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv), Z(CropRespv)) :
               OLP ?
    ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, Z(WindZ), Z(LeeZ), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                   Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv)) :
               N43P ?
    ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                   Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv)) :
    ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                   Z(VentSPv), Z(HumSPv), Z(Blackv))

FIXN = ["nu4", "U8_min_leak", "heat_ctrl", "hc_Kp", "hc_Ki", "vent_screen", "stomata_fix", "vapour_fix",
        "cond_heat", "f_cond", "vent_ctrl", "vc_hband", "pipe_lag", "black_scr"]
FIXV = [4.46e-5, 0.0, 1.0, 10.0, 10.0 / 1800, 1.0, 1.0, 1.0,
        1.0, 0.0, 1.0, 1.0, 0.0, 0.0]
if N43P     # Note 43 structure (all default-off engine flags, switched on here)
    append!(FIXN, ["hc_log_pipe", "vc_min_log", "soil_dyn"]); append!(FIXV, [1.0, 1.0, 1.0])
end
if N46P   # Notes 46-48: closed floor balance, solar split, own floor-deep conductance, internal mass;
                     # sd_tau fixed at 14 d (insensitive in SF, Note 47)
    append!(FIXN, ["floor_full", "solar_split", "deep_k", "mass_node", "sd_tau"])
    append!(FIXV, [1.0, 1.0, 1.0, 1.0, 14 * 86400.0])
end
if PRIORS == "n54"   # Note 54b: full lamp electrical power (lamp_Pel 117 W m-2 fixed from the AGC ReadMe)
    append!(FIXN, ["lamp_full"]); append!(FIXV, [1.0])
end
if OLP   # Note 58d: open loop -- logged lee/wind openings drive vent_split; lamp_full with screen loss
    FIXV[findfirst(==("vent_ctrl"), FIXN)] = 0.0
    append!(FIXN, ["lamp_full", "vent_split", "vent_exp", "vs_clee"]); append!(FIXV, [1.0, 1.0, 1.0, 1.0])
end
if N79   # Note 79: physical description fixed a priori (AGC documents, Notes 67-68, 78); bias terms removed
    append!(FIXN, ["phi2", "alpha5", "nu5", "gamma1", "n_pipes", "tau1", "psi2", "cover_area", "cover_cond", "co2_resp",
                   "hc_dT", "lamp_scr_loss", "s_rad", "dew_bl", "leak_shared", "rd_q10"])
    append!(FIXV, [5.4, 1005.0, 1000.0, 1.0, 1.0, 0.65 / 0.9, 35556.0, 1.6, 1.0, 1.0,
                   0.0, 0.0, 0.0, 1.0, 1.0, 1.0])     # dew_bl: Note 85; leak_shared, rd_q10: Note 86
    #  phi2 5.4 m mean air height | alpha5 c_p 1005 | nu5 = 0.10 x alpha6 vent area | rail pipe 1.0 m m-2 as ONE circuit
    #  (gamma1 1.0, n_pipes 1; PIPE2 folds the grow pipe in) | tau1 = 0.65/(1-eta1) -> PAR/global transmission 0.65 |
    #  psi2 = 12.8 g m-2 h-1 x alpha6 | cover 1.6 m2 per m2 floor | cover condensation and crop respiration ON
    allunique(FIXN) || error("duplicate FIXN entries")
end
const params_cl = update_params(params, FIXN, FIXV)

# ── free parameters ──────────────────────────────────────────────────────────
const CQ_MODEL = N79 ?
    ["gamma4", "vent_scr_perm", "k_ground", "nu1", "nu2", "C_soil", "k_deep", "C_mass", "h_mass", "f_mass",
     "lamp_Pel", "lamp_fcan", "vs_gwind", "vs_cwind", "vs_exp", "nu4", "h_cover_out", "air_cap_mult"] :
                 OLP ?
    ["gamma4", "vent_scr_perm", "hc_dT", "k_ground", "gamma1", "nu1", "nu2", "C_soil", "k_deep", "C_mass", "h_mass", "f_mass",
     "lamp_Pel", "lamp_scr_loss", "lamp_fcan", "vs_gwind", "vs_cwind", "vs_exp",
     (PRIORS == "n61" ? ["air_cap_mult", "s_rad"] : PRIORS == "n61b" ? ["air_cap_mult"] : String[])...] :
    vcat(["gamma4", "vent_exp", "vent_scr_perm", "hc_dT", "k_ground", "vc_gain", "vc_e0", "vc_tau", "vc_hum"],
                      PRIORS == "n37" ? ["vc_hband"] : PRIORS == "n39" ? ["vc_hband", "vc_min_night"] :
                      PRIORS == "n43" ? ["vc_hband", "gamma1", "sd_tau"] :
                      PRIORS == "n46" ? ["vc_hband", "gamma1", "nu1", "C_soil", "k_deep", "C_mass", "h_mass", "f_mass"] :
                      PRIORS == "n49b" ? ["vc_hband", "gamma1", "nu1", "nu2", "C_soil", "k_deep", "C_mass", "h_mass", "f_mass"] :
                      PRIORS == "n54" ? ["vc_hband", "gamma1", "nu1", "nu2", "C_soil", "k_deep", "C_mass", "h_mass", "f_mass", "lamp_fcan"] :
                      PRIORS == "n49a" ? ["vc_hband", "gamma1", "nu1", "nu2", "C_soil", "k_deep", "C_mass", "h_mass", "f_mass", "dday"] : String[])
const CO2N     = N79 && get(ENV, "CO2NIGHT", "1") == "1"   # Note 86: dark-hour CO2 as a third observable (n79 default on)
const CQ       = CO2N ? vcat(CQ_MODEL, ["sigma_T", "sigma_RH", "sigma_C"]) :
                 OLP ? vcat(CQ_MODEL, ["sigma_T", "sigma_RH"]) : vcat(CQ_MODEL, ["sigma_T", "sigma_RH", "sigma_U"])   # n58: no vent term
const NM       = length(CQ_MODEL)
const CQ_dict = Dict(                     # [min, centre, max]
    "gamma4"        => [20.0,  82.0,  300.0],
    "vent_exp"      => [0.3,   0.8,   1.5],
    "vent_scr_perm" => [0.0,   0.4,   1.0],
    "hc_dT"         => [-3.0,  0.0,   3.0],
    "k_ground"      => [1.0,   10.0,  30.0],
    "vc_gain"       => [0.03,  0.12,  0.5],     # 1/K   (Note 34: P-band ~9 K -> 0.11)
    "vc_e0"         => [-4.0, -2.0,   1.0],     # K     (Note 34: opening starts ~-2.5 K)
    "vc_tau"        => [120.0, 900.0, 7200.0],  # s     (Note 34: reaction <= 5 min, but hourly data)
    "vc_hum"        => [0.0,   0.12,  0.5],     # humidity-venting opening (Note 35: logged 0.15-0.22)
    "sigma_T"       => [0.2,   1.5,   6.0],
    "sigma_RH"      => [0.5,   6.0,   25.0],
    "sigma_U"       => [0.01,  0.08,  0.5],
    "sigma_C"       => [10.0,  60.0,  400.0])   # mg m-3, dark-hour CO2 (Note 86)
if PRIORS in ("n37", "n39", "n43", "n46", "n49a", "n49b", "n54")     # Note 37: soft humidity band removes the limit cycle; Note 34 controller identification
    CQ_dict["vc_gain"]  = [0.05,  0.12,  0.3]      # 1/K
    CQ_dict["vc_e0"]    = [-4.0, -2.0,   0.0]      # K
    CQ_dict["vc_tau"]   = [60.0,  600.0, 1800.0]   # s  (capped at 30 min)
    CQ_dict["vc_hum"]   = [0.0,   0.3,   1.0]
    CQ_dict["vc_hband"] = [0.5,   3.0,   6.0]      # g/m3 (was fixed at 1)
end
if PRIORS == "n39"     # Note 40: night-ventilation deficit; minimum night opening (engine flag, Note 38)
    CQ_dict["vc_min_night"] = [0.0, 0.03, 0.15]
end
if PRIORS == "n43"     # Note 43: pipe heat transfer (pipe length), deep-soil time constant, wider floor coupling
    CQ_dict["gamma1"]   = [0.6, 2.5, 6.0]                      # m/m2 (JSON 1.25; forward test favoured x2)
    CQ_dict["sd_tau"]   = [2 * 86400.0, 7 * 86400.0, 30 * 86400.0]   # s
    CQ_dict["k_ground"] = [1.0, 20.0, 80.0]                    # W/m2/K (was bound-pinned at 30)
end
if N46P     # physical priors (Notes 46-48): literature centres, physically sized ranges
    CQ_dict["gamma4"]   = [40.0,  82.0,  200.0]    # s/m   Vanthoor minimum stomatal resistance 82
    CQ_dict["nu1"]      = [0.4,   0.65,  1.0]      # -     vent discharge coefficient, Vanthoor 0.65
    CQ_dict["k_ground"] = [1.0,   8.0,   40.0]     # W/m2/K air <-> floor
    CQ_dict["C_soil"]   = [1.5e5, 6e5,   2e6]      # J/m2/K active floor layer ~0.1-1.3 m
    CQ_dict["k_deep"]   = [0.5,   2.0,   10.0]     # W/m2/K floor -> deep soil
    CQ_dict["gamma1"]   = [0.6,   2.5,   6.0]      # m/m2  pipe length (JSON 1.25; screens favour x2)
    CQ_dict["C_mass"]   = [1e4,   5e4,   3e5]      # J/m2/K internal mass (structure, slabs, water)
    CQ_dict["h_mass"]   = [1.0,   5.0,   20.0]     # W/m2/K air <-> internal mass
    CQ_dict["f_mass"]   = [0.0,   0.15,  0.35]     # -     share of outside radiation absorbed by the mass
end
if PRIORS in ("n49a", "n49b", "n54", "n58", "n61", "n61b")   # Note 50: physical caps; vent exchange per opening allowed lower; wind term free
    CQ_dict["nu1"]    = [0.3,  0.65, 1.0]
    CQ_dict["nu2"]    = [0.03, 0.10, 0.2]      # wind pressure coefficient of the roof-vent flow (JSON 0.1)
    CQ_dict["k_deep"] = [0.5,  1.5,  3.0]      # no large deep store (Note 49: ~5 m at k_deep 9)
    CQ_dict["f_mass"] = [0.0,  0.1,  0.2]      # inside solar absorption <= ~0.85 of outside
    CQ_dict["dday"]   = [0.0,  0.5,  2.0]      # K, daytime observation offset (n49a only; measurement model)
end
if OLP   # Note 58d: open-loop priors (lamp heat retained by day ~45-55, at night ~85-100 W m-2; windward vents)
    CQ_dict["nu2"]           = [0.03, 0.10, 0.5]
    CQ_dict["lamp_Pel"]      = [60.0, 117.0, 130.0]
    CQ_dict["lamp_scr_loss"] = [0.0,  0.4,  0.8]
    CQ_dict["lamp_fcan"]     = [0.05, 0.3,  0.9]
    CQ_dict["vs_gwind"]      = [0.5,  1.5,  4.0]
    CQ_dict["vs_cwind"]      = [0.5,  1.0,  10.0]
    CQ_dict["vs_exp"]        = [0.5,  0.9,  1.2]
end
if N61P   # Note 61: effective air heat capacity; sun-proportional term on the COMPARED T (measurement model)
    CQ_dict["air_cap_mult"] = [1.0, 4.0, 12.0]
    CQ_dict["s_rad"]        = [0.0, 0.1, 0.4]    # K per 100 W m-2 of outside global radiation
end
if PRIORS == "n61b"   # Note 62: no sensor term (s_rad stays 0 from the JSON); more absorbed sun may reach the air
    CQ_dict["f_mass"]       = [0.0, 0.25, 0.5]
    CQ_dict["air_cap_mult"] = [1.0, 6.0, 15.0]
    CQ_dict["h_mass"]       = [1.0, 5.0, 40.0]     # was capped at 20 in every fit since n46
end
if PRIORS == "n54"   # Note 54b: nu2 was pinned at 0.03 in n49b and OL is under-ventilated by day -> wider; lamp split free
    CQ_dict["nu2"]       = [0.03, 0.10, 0.5]
    CQ_dict["lamp_fcan"] = [0.1,  0.45, 0.7]   # share of lamp power absorbed by a closed canopy
end
if N79   # Note 78 prior table (approximate sub-models only)
    for (k, v) in ("nu1" => [0.3, 0.65, 1.0], "nu2" => [0.005, 0.05, 0.3], "vs_exp" => [0.5, 0.9, 1.2],
                   "vs_gwind" => [0.5, 1.5, 4.0], "vs_cwind" => [0.5, 1.0, 10.0], "nu4" => [1e-5, 4.5e-5, 1e-4],
                   "gamma4" => [40.0, 82.0, 200.0], "lamp_Pel" => [100.0, 110.0, 120.0], "lamp_fcan" => [0.05, 0.3, 0.9],
                   "k_ground" => [1.0, 8.0, 40.0], "C_soil" => [1.5e5, 6e5, 2e6], "k_deep" => [0.5, 1.5, 3.0],
                   "C_mass" => [2e4, 4e4, 8e4], "h_mass" => [1.0, 5.0, 40.0], "f_mass" => [0.0, 0.15, 0.5],
                   "h_cover_out" => [10.0, 25.0, 50.0], "air_cap_mult" => [1.0, 3.0, 12.0], "vent_scr_perm" => [0.0, 0.4, 1.0])
        CQ_dict[k] = v
    end
end
const ADDITIVE = Set(["hc_dT", "vc_e0"])        # jittered additively (centre 0 / negative)
CSupp(z) = all(CQ_dict[CQ[k]][1] < z[k] < CQ_dict[CQ[k]][3] for k in eachindex(z))

# ── solver with explicit tolerance (mirrors GreenhouseSim._run) ─────────────
function simdays(p, ds, u0x, rtol)
    us = Matrix{Float64}(undef, length(u0x), 1 + length(ds) * NPD); us[:, 1] = u0x
    u = copy(u0x); col = 2
    for d in ds
        ctx = GreenhouseSim.SimContext(p, weather_c, ctrl_x, t0 + 86400.0 * d, 86400.0 * d)
        sol = try
            GreenhouseSim.solve(GreenhouseSim.ODEProblem(GreenhouseSim.rhs!, u, tspan, ctx), GreenhouseSim.Tsit5();
                                saveat = teval, reltol = rtol)
        catch
            return nothing
        end
        length(sol.t) == NPD || return nothing
        S = Array(sol); us[:, col:col+NPD-1] = S; col += NPD; u = S[:, end]
    end
    return us
end

# ── blocks: 1 spin-up day + scored days, measured initial state, 8 states ───
struct Block
    days::Vector{Float64}; u0::Vector{Float64}
    binsT::Vector{Vector{Int}}; obsT::Vector{Float64}
    binsR::Vector{Vector{Int}}; obsR::Vector{Float64}
    binsU::Vector{Vector{Int}}; obsU::Vector{Float64}
    dmask::Vector{Float64}        # 1 where outside radiation > 50 W m-2 (for the n49a daytime observation offset)
    Ib::Vector{Float64}           # outside global radiation at the block columns (n61: s_rad measurement term)
    binsC::Vector{Vector{Int}}; obsC::Vector{Float64}   # Note 86: dark-hour CO2 bins (sun < 1, lamps off, valve closed)
end
const AT = A(Tair_obs); const AR = A(RH_obs); const AC = A(CO2_obs); const AU = A(U8log)
function make_block(a, b)
    d0 = max(a - 1, 0); ds = collect(Float64, d0:b)
    tstart = t0 + 86400.0 * d0
    Ta = AT(tstart); Ta = isnan(Ta) ? Tair_obs[findfirst(!isnan, Tair_obs)] : Ta
    Rh = AR(tstart); Rh = isnan(Rh) ? 80.0 : Rh
    Co = AC(tstart); Co = isnan(Co) ? 575.0 : Co
    Us = AU(tstart); Us = isnan(Us) ? 0.0 : Us
    u0b = [Ta, Ta, Co, Rh / 100 * Pws(Ta), 0.0, Ta, Us, Ta]
    if N43P     # deep-soil state starts at the mean measured air T of the (unscored) spin-up day
        ts_ = filter(!isnan, AT.(tstart .+ (0.0:300.0:86400.0))); push!(u0b, isempty(ts_) ? Ta : mean(ts_))
    end
    N46P && push!(u0b, u0b[end])      # internal mass starts at the same spin-up-day mean
    times = vcat(tstart, [t0 + d * 86400.0 + tj for d in ds for tj in teval])
    To = AT.(times); Ro = AR.(times); Uo = AU.(times); Co_ = AC.(times)
    darkC = (A(Ig_v).(times) .< 1.0) .& (A(U12v).(times) .< 0.01) .& (A(U10v).(times) .<= 0.0)
    bT = Vector{Int}[]; oT = Float64[]; bR = Vector{Int}[]; oR = Float64[]; bU = Vector{Int}[]; oU = Float64[]
    bC = Vector{Int}[]; oC = Float64[]
    for (k, d) in enumerate(ds)
        Int(d) < a && continue                                  # spin-up day not scored
        for hb in 0:23
            cols = [2 + (k - 1) * NPD + j for j in 0:NPD-1 if min(floor(Int, j * 300 / 3600), 23) == hb]
            okT = [c for c in cols if !isnan(To[c])]; okR = [c for c in cols if !isnan(Ro[c])]
            okU = [c for c in cols if !isnan(Uo[c])]
            length(okT) >= length(cols) / 2 && (push!(bT, okT); push!(oT, mean(To[okT])))
            length(okR) >= length(cols) / 2 && (push!(bR, okR); push!(oR, mean(Ro[okR])))
            length(okU) >= length(cols) / 2 && (push!(bU, okU); push!(oU, mean(Uo[okU])))
            okC = [c for c in cols if !isnan(Co_[c]) && darkC[c]]                  # whole hour dark, lamps off, valve shut
            length(okC) == length(cols) && (push!(bC, okC); push!(oC, mean(Co_[okC])))
        end
    end
    return Block(ds, u0b, bT, oT, bR, oR, bU, oU, Float64.(A(Ig_v).(times) .> 50.0), A(Ig_v).(times), bC, oC)
end
const TB  = [make_block(a, b) for (a, b) in TRAIN_BLOCKS]
const NBT = sum(length(b.obsT) for b in TB); const NBR = sum(length(b.obsR) for b in TB)
const NBU = sum(length(b.obsU) for b in TB)
const NBC = sum(length(b.obsC) for b in TB)

function block_sse(p, b::Block, rtol)
    us = simdays(p, b.days, b.u0, rtol); us === nothing && return (Inf, Inf, Inf, Inf)
    T2 = us[2, :] .+ (hasproperty(p, :dday) ? p.dday : 0.0) .* b.dmask   # n49a: compared T = air T + daytime offset
    hasproperty(p, :s_rad) && p.s_rad != 0.0 && (T2 = T2 .+ p.s_rad .* b.Ib ./ 100)   # n61: sun-heated sensor term
    RH = rhf.(us[2, :], us[4, :]); U8 = clamp.(us[7, :], 0.0, 1.0)
    sT = 0.0; sR = 0.0; sU = 0.0; sC = 0.0
    if !isempty(b.binsC); C1s = us[3, :]; for k in eachindex(b.binsC); sC += (binmean(C1s, b.binsC[k]) - b.obsC[k])^2; end; end
    for k in eachindex(b.binsT); sT += (binmean(T2, b.binsT[k]) - b.obsT[k])^2; end
    for k in eachindex(b.binsR); sR += (binmean(RH, b.binsR[k]) - b.obsR[k])^2; end
    for k in eachindex(b.binsU); sU += (binmean(U8, b.binsU[k]) - b.obsU[k])^2; end
    f(s) = isfinite(s) ? s : Inf
    return (f(sT), f(sR), f(sU), f(sC))
end
function sse_all(x, rtol)
    p = update_params(params_cl, CQ_MODEL, x)
    sT = 0.0; sR = 0.0; sU = 0.0; sC = 0.0
    for b in TB
        a, c, e, g = block_sse(p, b, rtol); (isfinite(a) && isfinite(c) && isfinite(e) && isfinite(g)) || return (Inf, Inf, Inf, Inf)
        sT += a; sR += c; sU += e; sC += g
    end
    return (sT, sR, sU, sC)
end
function cenergy(z; rtol = RELTOL)
    σT = z[NM+1]; σR = z[NM+2]
    sT, sR, sU, sC = sse_all(z[1:NM], rtol); isfinite(sT) || return Inf
    CO2N && return NBT * log(σT) + sT / (2σT^2) + NBR * log(σR) + sR / (2σR^2) + NBC * log(z[NM+3]) + sC / (2z[NM+3]^2)
    OLP && return NBT * log(σT) + sT / (2σT^2) + NBR * log(σR) + sR / (2σR^2)   # n58: vents are inputs
    σU = z[NM+3]
    return NBT * log(σT) + sT / (2σT^2) + NBR * log(σR) + sR / (2σR^2) + NBU * log(σU) + sU / (2σU^2)
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
        f["chain"] = CP.obj.Output[1:rows, :]; f["QoI"] = CQ; f["design"] = DESIGN; f["seed"] = CSEED; f["reltol"] = RELTOL
    end
    mv(CCKPT * ".tmp", CCKPT; force = true)
end
function csupp(z)
    CP.calls += 1; creport()
    if time() - CP.tc >= 1800; CP.tc = time(); ccheckpoint(); end
    return CSupp(z)
end
cenergy_t(z) = (U = cenergy(z); U < CP.Umin && (CP.Umin = U); U)

centre = [CQ_dict[q][2] for q in CQ]   # (overrides the joint script's `centre`)
function cjitter(c)
    for _ in 1:1000
        x = [CQ[k] in ADDITIVE ? c[k] + 0.2 * randn() : c[k] * (1 + 0.05 * randn()) for k in eachindex(c)]
        CSupp(x) && return x
    end
    error("no start point in support")
end

println("\n", "="^78)
@printf("CLOSED-LOOP climate calibration -- design=%s  priors=%s  mode=%s  seed=%d  reltol=%.0e\n", DESIGN, PRIORS, CMODE, CSEED, RELTOL)
println("  TRAIN days: ", TRAIN_BLOCKS, "   VAL: ", VAL_BLOCKS)
@printf("  TRAIN hourly bins: T %d, RH %d, VENT %d, dark CO2 %d%s  (%d simulated days per step incl. spin-up)\n",
        NBT, NBR, NBU, NBC, CO2N ? " (in likelihood)" : " (not used)", sum(length(b.days) for b in TB))
println("="^78)

# ── evaluation on the full season (closed loop from day 0) ──────────────────
const u0_8 = vcat(u0, U8log[findfirst(!isnan, U8log)], u0[2])
const u0_cl = N43P ? vcat(u0_8, mean(filter(!isnan, AT.(t0 .+ (0.0:300.0:86400.0))))) : u0_8
N46P && push!(u0_cl, u0_cl[end])      # internal mass (Note 48)
function full_eval(x; rtol = RELTOL)
    us = simdays(update_params(params_cl, CQ_MODEL, x), days, u0_cl, rtol)
    us === nothing && error("full-season solve failed")
    dd = "dday" in CQ_MODEL ? x[findfirst(==("dday"), CQ_MODEL)] : 0.0       # n49a: compared T includes the offset
    sr = "s_rad" in CQ_MODEL ? x[findfirst(==("s_rad"), CQ_MODEL)] : 0.0     # n61: compared T includes s_rad * I/100
    Uout = OLP ? AU.(predtimes) : clamp.(us[7, :], 0.0, 1.0)        # n58: vents are logged inputs (report the log)
    return us[2, :] .+ dd .* (A(Ig_v).(predtimes) .> 50.0) .+ sr .* A(Ig_v).(predtimes) ./ 100, rhf.(us[2, :], us[4, :]), Uout
end
cd_  = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
nite = A(Ig_v).(predtimes) .< 1.0
dayt = A(Ig_v).(predtimes) .> 50.0
U8obs = AU.(predtimes)
function hourly_diff(Um, mask)          # hourly-mean (model - logged) over masked columns, >= 6 samples/hour
    acc = Dict{Tuple{Int,Int},Tuple{Float64,Float64,Int}}()
    for c in 2:NCOL
        mask[c] || continue
        k = (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600))
        a, b, m = get(acc, k, (0.0, 0.0, 0)); acc[k] = (a + Um[c], b + U8obs[c], m + 1)
    end
    [a / m - b / m for (a, b, m) in values(acc) if m >= 6]
end
function eval_table(T2, RH, U8m; label = "")
    okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g); okU = .!isnan.(U8obs); okU[1] = false
    VPo = RHobs_g ./ 100 .* Pws.(Tobs_g) ./ 100; VPm = RH ./ 100 .* Pws.(T2) ./ 100
    sets = [("TRAIN", TRAIN_DAYS), ("VAL", VAL_DAYS), ("TEST", TEST_DAYS), ("COMMON TEST", COMMON_TEST)]
    @printf("\n  %s\n  %-12s %6s %7s %7s %7s %7s %7s | %7s %7s %7s %7s\n", label, "set", "RMSE", "bias", "night", "day",
            "VPday", "VPnite", "U8 RMSE", "U8 bias", "U8 day", "U8 nite")
    for (nm_, ds) in sets
        ins = in.(cd_, Ref(Set(ds))); m = okT .& ins; mv_ = okV .& ins
        mb(v, o, k) = any(k) ? mean(v[k] .- o[k]) : NaN
        hd = hourly_diff(U8m, okU .& ins); hdd = hourly_diff(U8m, okU .& ins .& dayt); hdn = hourly_diff(U8m, okU .& ins .& nite)
        mm(v) = isempty(v) ? NaN : mean(v)
        @printf("  %-12s %6.2f %+7.2f %+7.2f %+7.2f %+7.2f %+7.2f | %7.3f %+7.3f %+7.3f %+7.3f\n", nm_,
                sqrt(mean((T2[m] .- Tobs_g[m]) .^ 2)), mb(T2, Tobs_g, m), mb(T2, Tobs_g, m .& nite), mb(T2, Tobs_g, m .& dayt),
                mb(VPm, VPo, mv_ .& dayt), mb(VPm, VPo, mv_ .& nite),
                isempty(hd) ? NaN : sqrt(mean(abs2, hd)), mm(hd), mm(hdd), mm(hdn))
    end
    println("  U8 columns: hourly-mean vent opening, model minus logged (bias by day / night).")
    # within-hour variability (Note 37 limit-cycle indicator), COMMON TEST, model vs observed
    function wsd_(v, mask)
        acc = Dict{Tuple{Int,Int},Vector{Float64}}()
        for c in 2:NCOL
            mask[c] || continue
            push!(get!(acc, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Float64[]), v[c])
        end
        s = [std(x) for x in values(acc) if length(x) >= 9]; isempty(s) ? NaN : mean(s)
    end
    com = in.(cd_, Ref(Set(COMMON_TEST)))
    mT = okT .& com; mU = okU .& com
    @printf("  within-hour SD, COMMON TEST (model / observed):  T day %.3f / %.3f   T night %.3f / %.3f   U8 day %.3f / %.3f   U8 night %.3f / %.3f\n",
            wsd_(T2, mT .& dayt), wsd_(Tobs_g, mT .& dayt), wsd_(T2, mT .& nite), wsd_(Tobs_g, mT .& nite),
            wsd_(U8m, mU .& dayt), wsd_(U8obs, mU .& dayt), wsd_(U8m, mU .& nite), wsd_(U8obs, mU .& nite))
end

if CMODE == "smoke"
    # [smoke 1] tolerance check: energy level AND energy differences between two points
    xa = copy(centre); Random.seed!(CSEED + 100); xb = cjitter(centre); Random.seed!(CSEED)
    println("\n[smoke 1] tolerance check (energy at prior centre xa and at a 5%-jittered point xb)")
    cenergy(xa; rtol = 1e-3)                                     # compile
    @printf("  %-8s %12s %12s %12s %10s\n", "reltol", "U(xa)", "U(xb)", "U(xb)-U(xa)", "ms/eval")
    res = Dict{Float64,Tuple{Float64,Float64}}()
    for rt in (1e-3, 1e-4, 1e-5)
        te = @elapsed (Ua = cenergy(xa; rtol = rt)); Ub = cenergy(xb; rtol = rt); res[rt] = (Ua, Ub)
        @printf("  %-8.0e %12.2f %12.2f %12.2f %10.0f\n", rt, Ua, Ub, Ub - Ua, 1e3te)
    end
    d3 = (res[1e-3][2] - res[1e-3][1]) - (res[1e-5][2] - res[1e-5][1])
    d4 = (res[1e-4][2] - res[1e-4][1]) - (res[1e-5][2] - res[1e-5][1])
    @printf("  error in the energy DIFFERENCE vs 1e-5:  1e-3 -> %+.2f   1e-4 -> %+.2f\n", d3, d4)
    @printf("  verdict: %s\n", abs(d3) < 1.0 ? "reltol 1e-3 is adequate (|error| < 1 energy unit)" :
            abs(d4) < 1.0 ? "use RELTOL=1e-4 for the full run" : "use RELTOL=1e-5 for the full run")
    isfinite(res[1e-3][1]) || error("energy not finite at prior centre")
    T3, _, U3 = full_eval(xa[1:NM]; rtol = 1e-3); T5, _, U5 = full_eval(xa[1:NM]; rtol = 1e-5)
    @printf("  full season at centre, 1e-3 vs 1e-5: max |ΔT| %.2f K, mean |ΔT| %.3f K, max |ΔU8| %.3f\n",
            maximum(abs.(T3 .- T5)), mean(abs.(T3 .- T5)), maximum(abs.(U3 .- U5)))

    # [smoke 1b] energy split at the centre (which term dominates?)
    sT, sR, sU, sC = sse_all(xa[1:NM], RELTOL)
    @printf("[smoke 1b] at centre: RMS hourly misfit T %.2f K, RH %.2f %%, vent %.3f, dark CO2 %.0f mg/m3  (bins %d / %d / %d / %d)\n",
            sqrt(sT / NBT), sqrt(sR / NBR), sqrt(sU / NBU), NBC > 0 ? sqrt(sC / NBC) : NaN, NBT, NBR, NBU, NBC)

    # [smoke 2] short chain + projected cost
    obj = jtwalk(n = length(CQ), U = cenergy_t, Supp = csupp)
    CP.total = CSMOKE; CP.calls = 0; CP.t0 = CP.tl = CP.tc = time(); CP.obj = obj
    tw = @elapsed Run!(obj, T = CSMOKE, x0 = cjitter(centre), xp0 = cjitter(centre))
    creport(true); println()
    O = obj.Output; acc = mean(any(O[k+1, 1:end-1] .!= O[k, 1:end-1]) for k in 1:size(O, 1)-1)
    @printf("[smoke 2] %d steps in %.1f s (%.0f ms/step at reltol %.0e), moved-fraction %.2f, U %.1f -> %.1f\n",
            CSMOKE, tw, 1e3tw / CSMOKE, RELTOL, acc, O[1, end], O[end, end])
    @printf("          projected full run (%d steps): %.1f h at reltol %.0e\n", CCHAIN, tw / CSMOKE * CCHAIN / 3600, RELTOL)

    # [smoke 3] full-season check at the prior centre
    T2c, RHc, U8c = full_eval(centre[1:NM]); eval_table(T2c, RHc, U8c; label = "Full-season check at the PRIOR CENTRE (not calibrated):")
    if N79   # [smoke 4] the two new processes at the prior centre: CO2 and night vapour (checks, not in the likelihood)
        Cobs = A(CO2_obs).(predtimes); RHo = RHobs_g
        println("\n[smoke 4] new processes at the prior centre, full season (CO2 and RH are checks / likelihood as before)")
        @printf("  %-28s %9s %9s %9s %9s %9s\n", "variant", "CO2 RMSE", "CO2 bias", "CO2 nite", "RH nite", "RH day")
        for (lab, fx) in (("both on", Pair{String,Float64}[]), ("crop respiration off", ["co2_resp" => 0.0]),
                          ("cover condensation off", ["cover_cond" => 0.0]), ("both off", ["co2_resp" => 0.0, "cover_cond" => 0.0]))
            pv = update_params(params_cl, CQ_MODEL, centre[1:NM])
            isempty(fx) || (pv = update_params(pv, String[first(q) for q in fx], Float64[last(q) for q in fx]))
            us = simdays(pv, days, u0_cl, RELTOL); us === nothing && (println("  $lab: solve failed"); continue)
            C = us[3, :]; RHm = rhf.(us[2, :], us[4, :])
            okc = .!isnan.(Cobs); okc[1] = false; okr = .!isnan.(RHo); okr[1] = false
            mbx(v, o, k) = any(k) ? mean(v[k] .- o[k]) : NaN
            @printf("  %-28s %9.0f %+9.0f %+9.0f %+9.2f %+9.2f\n", lab, sqrt(mean((C[okc] .- Cobs[okc]) .^ 2)), mbx(C, Cobs, okc),
                    mbx(C, Cobs, okc .& nite), mbx(RHm, RHo, okr .& nite), mbx(RHm, RHo, okr .& dayt))
        end
        println("  CO2 in mg m-3 (x ~0.55 -> ppm); RH in % points; bias = model - observed.")
    end
    println("\nSmoke done. Paste this output before launching `full`.")
    exit(0)
end

# ── full run ─────────────────────────────────────────────────────────────────
if CMODE == "full"
obj = jtwalk(n = length(CQ), U = cenergy_t, Supp = csupp)
x0c, xp0c = cjitter(centre), cjitter(centre)
CP.total = CCHAIN; CP.calls = 0; CP.t0 = CP.tl = CP.tc = time(); CP.obj = obj
@printf("Progress every %.0f s; checkpoint every 30 min -> %s\n", CPSEC, CCKPT)
@time Run!(obj, T = CCHAIN, x0 = x0c, xp0 = xp0c)
creport(true); println()

post = obj.Output[CBURN+1:end, 1:end-1]; pm = vec(mean(post, dims = 1))
T2m, RHm, U8m = full_eval(pm[1:NM])
JLD2.@save joinpath(POSTDIR, CTAG * ".jld2") chain=obj.Output QoI=CQ burn_in=CBURN design=DESIGN seed=CSEED reltol=RELTOL train=TRAIN_BLOCKS val=VAL_BLOCKS fixn=FIXN fixv=FIXV priors=PRIORS T2m RHm U8m
println("Chain saved: ", joinpath(POSTDIR, CTAG * ".jld2"))

println("\n--- posterior ---")
for (i, q) in enumerate(CQ)
    lo, _, hi = CQ_dict[q]; qs = quantile(post[:, i], [0.025, 0.975]); span = hi - lo
    flag = qs[1] - lo < 0.02span ? "AT LOWER BOUND" : hi - qs[2] < 0.02span ? "AT UPPER BOUND" : ""
    @printf("  %-14s mean %9.4g   95%% [%9.4g, %9.4g]  %s\n", q, pm[i], qs[1], qs[2], flag)
end
acc = mean(any(obj.Output[k+1, 1:end-1] .!= obj.Output[k, 1:end-1]) for k in 1:size(obj.Output, 1)-1)
@printf("moved-fraction %.3f\n", acc)
eval_table(T2m, RHm, U8m; label = "Full-season evaluation at the POSTERIOR MEAN (closed-loop vents):")
end   # CMODE == "full"
