# =============================================================================
# equivalence_n79.jl -- does GreenhouseTwin (flag-free) reproduce V2.2 with the n79 flags?
# =============================================================================
# Run from invernaderosN79/GreenhouseTwin (a few minutes; no MCMC):
#     julia --project=. test/equivalence_n79.jl 2>&1 | tee equivalence_n79.log
# Optional: CROPTAG=n94 (default; which V2.2 crop series drives both climate models), NRAND=2000 (random rhs! points).
#
# Loads V2.2 through its own n79 set-up (test/run_climate_calibration_cl.jl in CL_MODE=setup, PRIORS=n79, SOLARC=1,
# CO2NIGHT=1, PIPE2=1) and GreenhouseTwin side by side, then compares:
#   1 inputs built independently by both codes (every array, bit for bit)
#   2 parameters (every GreenhouseTwin parameter against the V2.2 n79 value at the prior centre)
#   3 rhs! at random states and times (also exactly on the 5-min grid), relative to each derivative's scale
#   4 TRAIN blocks: bins, initial states, trajectories at reltol 1e-4
#   5 climate energy at the prior centre and at a jittered point
#   6 full season at reltol 1e-3 (both codes, timed)
#   7 solar budget closure
#   8 crop at the n105 prior centre: series and energy
#   9 rhs! timing and allocations
# Tolerances: inputs, parameters, bins exact; rhs! 1e-10 relative. Trajectories and energies are judged against the
# solver's own sensitivity: the model has discontinuous inputs and switches (PI controller, max()), so the adaptive
# solver amplifies rounding-level differences; V2.2 against itself with the initial state (or one parameter) scaled by
# 1 + 1e-14 gives the noise floor, and a check passes when |new - V2.2| <= 5 x floor. (The flag-free code books a few
# cancelling V2.2 terms once, e.g. lamp NIR added then removed, so derivatives agree to rounding, not bit for bit.)
# =============================================================================
const EQ_HERE = @__DIR__
const EQ_INV  = normpath(joinpath(EQ_HERE, "..", "..", ".."))
const EQ_OLD  = joinpath(EQ_INV, "invernaderosV2.2", "GreenhouseSim")
const EQ_CROPTAG = get(ENV, "CROPTAG", "n94")
const EQ_NRAND = parse(Int, get(ENV, "NRAND", "2000"))
for (k, v) in ("PRIORS" => "n79", "CL_MODE" => "setup", "CROPTAG" => EQ_CROPTAG, "SOLARC" => "1", "CO2NIGHT" => "1",
               "DESIGN" => "blocks", "PIPE2" => "1")
    ENV[k] = v
end
println("== loading V2.2 (n79 set-up script, CL_MODE=setup) ==")
include(joinpath(EQ_OLD, "test", "run_climate_calibration_cl.jl"))      # defines the V2.2 globals (Main)
println("== loading GreenhouseTwin ==")
include(joinpath(EQ_HERE, "..", "src", "GreenhouseTwin.jl"))
const GT = GreenhouseTwin
using Printf, Random, Statistics, BenchmarkTools

const EQ_NP = Ref(0); const EQ_NF = Ref(0)
eqcheck(name, ok) = (ok ? (EQ_NP[] += 1) : (EQ_NF[] += 1); @printf("%-92s %s\n", name, ok ? "PASS" : "FAIL"); ok)
maxabs(a, b) = maximum(abs.(filter(!isnan, a .- b)); init = 0.0)

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 1. inputs (both codes read the AGC files independently) ==")
const EQ_CS = joinpath(EQ_OLD, "test", "posteriors", "AiCU_crop_series_$(EQ_CROPTAG).csv")
const SN = GT.load_agc_season(EQ_INV; crop_series = EQ_CS)
for (nm, a, b) in (("time stamps", SN.ts, ts), ("global radiation", SN.Ig, Ig_v), ("sky T", SN.Tsky, Tsky_v),
                   ("outside T", SN.Tout, Tout_v), ("wind", SN.Wind, Wind_v), ("outside VP", SN.VPout, VPout_v),
                   ("energy screen", SN.U1, U1v), ("blackout screen", SN.Black, Blackv), ("lee vent", SN.Lee, LeeZ),
                   ("wind vent", SN.WindVent, WindZ), ("CO2 valve (dosage-scaled)", SN.U10, U10v),
                   ("heating setpoint", SN.U11, U11v), ("lamps", SN.U12, U12v), ("pipe T (n79 rule)", SN.Tpipe, Tpipe_meas),
                   ("LAI", SN.LAI, LAI_row), ("crop respiration", SN.CropResp, CropRespv), ("crop fresh mass", SN.CropFW, CropFWv),
                   ("obs air T", SN.Tair_obs, Tair_obs), ("obs RH", SN.RH_obs, RH_obs), ("obs CO2", SN.CO2_obs, CO2_obs))
    eqcheck(@sprintf("  %-28s identical (n = %d, max |diff| %.3g)", nm, length(a), length(a) == length(b) ? maxabs(a, b) : NaN),
            length(a) == length(b) && isequal(a, b))
end
eqcheck("  season start t0 identical", SN.t0 == t0)

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 2. parameters (GreenhouseTwin JSON vs V2.2 n79 at the prior centre) ==")
const CJ = joinpath(EQ_HERE, "..", "configfiles", "climate_n79.json")
const pN = GT.load_params(CJ)
const pO = update_params(params_cl, CQ_MODEL, centre[1:NM])
let bad = String[], newn = String[]
    for k in keys(pN)
        if hasproperty(pO, k)
            getproperty(pO, k) == getproperty(pN, k) || push!(bad, @sprintf("%s: V2.2 %.17g, new %.17g", k, getproperty(pO, k), getproperty(pN, k)))
        else
            push!(newn, String(k))
        end
    end
    foreach(x -> println("    differs: ", x), bad)
    println("    names only in GreenhouseTwin: ", join(newn, ", "))
    eqcheck(@sprintf("  %d shared parameters have identical values", length(keys(pN)) - length(newn)), isempty(bad))
    eqcheck("  new names are the constants that V2.2 hard-coded (k_lai_photo 0.7, g_s 0.3, g_m 0.25) + noise",
            sort(newn) == sort(["k_lai_photo", "g_s", "g_m", "sigma_T", "sigma_RH", "sigma_C"]) &&
            pN.k_lai_photo == 0.7 && pN.g_s == 0.3 && pN.g_m == 0.25)
    eqcheck("  free-parameter list and prior centres = V2.2 n79", GT.CLIMATE_QOI == CQ_MODEL &&
            [getproperty(pN, Symbol(q)) for q in CQ_MODEL] == centre[1:NM])
    fx = Dict(zip(FIXN, FIXV))
    eqcheck("  V2.2 n79 flags as expected (solar_cons, par_below, cover_cond, co2_resp, leak_stack, ... on)",
            all(fx[k] == 1.0 for k in ("solar_cons", "solar_mass_share", "par_below", "cover_cond", "co2_resp", "leak_stack",
                                         "leak_shared", "rd_q10", "dew_bl", "crop_cap", "lamp_full", "vent_split", "floor_full",
                                         "solar_split", "deep_k", "mass_node", "soil_dyn", "heat_ctrl", "hc_log_pipe",
                                         "vent_screen", "stomata_fix", "vapour_fix")) &&
            fx["vent_ctrl"] == 0.0 && fx["pipe_lag"] == 0.0 && fx["black_scr"] == 0.0 && fx["f_cond"] == 0.0 &&
            pO.use_cover_node == 1.0 && pO.phys_light == 1.0 && pO.photo_lai == 1.0)
end

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 3. rhs! at random states and times ==")
const IMAP = [1, 2, 3, 4, 5, 6, 9, 10]          # V2.2 state index of each GreenhouseTwin state
old_u(u) = [u[1], u[2], u[3], u[4], u[5], u[6], 0.0, u[2], u[7], u[8]]   # u7 (vent state), u8 (pipe lag) unused in n79
let rng = MersenneTwister(79), DO = zeros(EQ_NRAND, 8), DN = zeros(EQ_NRAND, 8)
    duo = zeros(10); dun = zeros(8)
    for k in 1:EQ_NRAND
        d = rand(rng, 0:114); t = isodd(k) ? 300.0 * rand(rng, 0:287) : 86400.0 * rand(rng)
        T2 = 285 + 20rand(rng); T1 = T2 + 3 * randn(rng); RH = 40 + 58rand(rng)
        u = [T1, T2, 400 + 1100rand(rng), RH / 100 * Pws(T2), 3000 * rand(rng) * (rand(rng) < 0.3 ? -1 : 1),
             T2 - 5 + 8rand(rng), 285 + 10rand(rng), T2 - 3 + 6rand(rng)]
        co = GreenhouseSim.SimContext(pO, weather_c, ctrl_x, t0 + 86400.0 * d, 86400.0 * d)
        cn = GT.SimContext(pN, SN.weather, SN.controls, SN.t0 + 86400.0 * d, 86400.0 * d)
        GreenhouseSim.rhs!(duo, old_u(u), co, t); GT.rhs!(dun, u, cn, t)
        DO[k, :] = duo[IMAP]; DN[k, :] = dun
    end
    for (j, nm) in enumerate(("dT1", "dT2", "dC1", "dV1", "de_int", "dT5", "dTd", "dTm"))
        sc = maximum(abs.(DO[:, j])); err = maximum(abs.(DO[:, j] .- DN[:, j])) / max(sc, 1e-300)
        eqcheck(@sprintf("  %-7s max |diff| / max |value| = %.2e  (scale %.3g)", nm, err, sc), err < 1e-10)
    end
end

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 4. TRAIN blocks ==")
const BN = GT.climate_blocks(SN)
eqcheck("  same number of blocks and observation bins (T, RH, dark CO2)",
        length(BN) == length(TB) && all(b.binsT == o.binsT && b.binsR == o.binsR && b.binsC == o.binsC for (b, o) in zip(BN, TB)))
eqcheck("  identical bin means (T, RH, CO2)", all(isequal(b.obsT, o.obsT) && isequal(b.obsR, o.obsR) && isequal(b.obsC, o.obsC)
                                               for (b, o) in zip(BN, TB)))
eqcheck("  identical initial states", all(b.u0 == o.u0[IMAP] for (b, o) in zip(BN, TB)))
let worst = zeros(4), floor_ = zeros(4)
    for (b, o) in zip(BN, TB)
        uo = simdays(pO, o.days, o.u0, 1e-4); un = GT.simdays(pN, SN, b.days, b.u0, 1e-4)
        up = simdays(pO, o.days, o.u0 .* (1 + 1e-14), 1e-4)                  # V2.2 against itself: solver noise floor
        for (j, r) in enumerate((2, 1, 3, 4))
            worst[j] = max(worst[j], maximum(abs.(uo[IMAP[r], :] .- un[r, :])))
            floor_[j] = max(floor_[j], maximum(abs.(uo[IMAP[r], :] .- up[IMAP[r], :])))
        end
    end
    for (j, nm) in enumerate(("T2 [K]", "T1 [K]", "C1 [mg m-3]", "V1 [Pa]"))
        eqcheck(@sprintf("  %-12s |new - V2.2| %.2e   noise floor |V2.2 - V2.2(u0 x (1+1e-14))| %.2e", nm, worst[j], floor_[j]),
                worst[j] <= 5 * floor_[j] + 1e-9)
    end
end

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 5. climate energy ==")
let zc = vcat([getproperty(pN, Symbol(q)) for q in GT.CLIMATE_QOI], pN.sigma_T, pN.sigma_RH, pN.sigma_C)
    eqcheck("  energy argument (prior centre) = V2.2 centre vector", zc == centre)
    # The likelihood is numerically noisy (discontinuous inputs and switches; see the header): compare ENSEMBLES of
    # energies over tiny perturbations (gamma4 x (1 + k 1e-14), k = 0..NENS-1) at reltol 1e-7, where the noise is small.
    # (At 1e-6 a systematic difference of ~0.5 appears at some points; it shrinks to ~0.03 at 1e-7: integration error
    # interacting with the arithmetic order, not a model difference.)
    nens = parse(Int, get(ENV, "NENS", "4"))
    for (lab, z0) in (("prior centre", centre), ("5%-jittered point", (Random.seed!(7); cjitter(centre))))
        zs = [(z = copy(z0); z[1] *= 1 + k * 1e-14; z) for k in 0:nens-1]
        Uo = [cenergy(z; rtol = 1e-7) for z in zs]; Un = [GT.climate_energy(z, pN, SN, BN; rtol = 1e-7) for z in zs]
        se = sqrt(var(Uo) / nens + var(Un) / nens)
        eqcheck(@sprintf("  U at the %s (reltol 1e-7, %d perturbations): V2.2 %.3f ± %.3f, new %.3f ± %.3f, diff %+.3f (3 se %.3f)",
                         lab, nens, mean(Uo), std(Uo), mean(Un), std(Un), mean(Un) - mean(Uo), 3se),
                abs(mean(Un) - mean(Uo)) <= 3se + 1e-6 * abs(mean(Uo)))
    end
    Uo4 = [cenergy((z = copy(centre); z[1] *= 1 + k * 1e-14; z); rtol = 1e-4) for k in 0:nens-1]
    @printf("    information: at reltol 1e-4 (calibration default) the V2.2 energy at the centre has numerical noise sd %.2f (range %.2f)\n",
            std(Uo4), -(-(extrema(Uo4)...)))
end

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 6. full season (115 days, reltol 1e-3) ==")
eqcheck("  day-0 initial state identical", SN.u0 == u0_cl[IMAP])
let
    simdays(pO, [0.0], u0_cl, 1e-3); GT.simdays(pN, SN, [0.0], SN.u0, 1e-3)      # compile
    to = @elapsed (uo = simdays(pO, days, u0_cl, 1e-3)); tn = @elapsed (un = GT.simdays(pN, SN, days, SN.u0, 1e-3))
    up = simdays(pO, days, u0_cl .* (1 + 1e-14), 1e-3)
    for (r, nm) in ((2, "T2 [K]"), (3, "C1 [mg m-3]"), (4, "V1 [Pa]"))
        e = maximum(abs.(uo[r, :] .- un[r, :])); b = maximum(abs.(uo[r, :] .- up[r, :]))
        eqcheck(@sprintf("  %-12s |new - V2.2| %.2e   noise floor %.2e", nm, e, b), e <= 5b + 1e-9)
    end
    @printf("    time: V2.2 %.1f s, new %.1f s\n", to, tn)
    ok = .!isnan.(Tobs_g); ok[1] = false
    @printf("    hourly-grid T RMSE vs observed at the prior centre: V2.2 %.4f K, new %.4f K\n",
            sqrt(mean((uo[2, ok] .- Tobs_g[ok]) .^ 2)), sqrt(mean((un[2, ok] .- Tobs_g[ok]) .^ 2)))
end

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 7. solar budget ==")
let rng = MersenneTwister(5), worst = 0.0
    for _ in 1:5000
        I = 1000rand(rng); s = GT.solar_column(I, rand(rng), rand(rng), 6rand(rng), pN)
        worst = max(worst, abs(s.Qcan + s.Qflr + s.Qair + s.Qcov + s.Qout - I))
    end
    eqcheck(@sprintf("  canopy + floor + air + cover + sky = Iglob (5000 random states, max error %.1e W m-2)", worst), worst < 1e-9)
    so = GreenhouseSim.solar_column(500.0, 0.7, 0.3, 2.3, pO); sn = GT.solar_column(500.0, 0.7, 0.3, 2.3, pN)
    eqcheck("  solar fluxes identical to V2.2 (one state, all fields)",
            all(getfield(sn, f) == getfield(so, f) for f in fieldnames(GT.SolarFluxes)))
end

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 8. crop (n105 set-up, prior centre) ==")
let
    # V2.2: as test/calibrate_canopy_growth.jl with PARBELOW=1
    D = joinpath(EQ_INV, "invernaderos2.0", "data")
    pco = load_params(joinpath(EQ_INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json"),
                      joinpath(EQ_INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json"))
    pco = update_params(pco, ["par_below", "tau1"], [1.0, 0.65 / 0.9])
    gpo = load_growth_params(joinpath(EQ_OLD, "configfiles", "cucumber_growth.json"))
    gpo = update_params(gpo, ["set_start_day", "harvest_W", "gen_clampT", "gen_parlight", "gen_beta", "gen_kappa", "gen_A",
                              "gen_V", "gen_zclamp", "stem_density", "gen_canopy", "gen_doy0", "gen_fL", "gen_flower_node"],
                        [15.0, 16.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.6, 1.0, 226.0, 0.59, 5.0])
    sun_s  = load_measured_climate_par(joinpath(D, "AiCU", "Greenhouse_climate.csv"), joinpath(D, "meteo.csv"); transmission = 1.0, lamp_umol = 0.0)
    lamp_s = load_measured_climate_par(joinpath(D, "AiCU", "Greenhouse_climate.csv"), joinpath(D, "meteo.csv"); transmission = 0.0, lamp_umol = 187.0)
    seaso  = MeasuredSeason(sun_s.ndays, [0.65 .* S .+ L for (S, L) in zip(sun_s.day_samples, lamp_s.day_samples)], sun_s.dt)
    for (W, S, L) in zip(seaso.day_samples, sun_s.day_samples, lamp_s.day_samples)
        W[1, :] = S[1, :]; W[2, :] = S[2, :]           # rows 1-2: T and CO2 (not scaled)
    end
    keepo  = load_fruit_keep(joinpath(D, "AiCU", "CropManagement.csv"); ndays = seaso.ndays, as_removal = true)
    V2 = CUCUMBER_PARAMS_V2
    genp(W_pot, E, b1, V20) = (AT = W_pot * V2.s * (V2.p - 1) / E;
        GenericParameters(s = V2.s, p = V2.p, peak = V2.peak, A0 = -10.0 * AT, AT = AT, kappa = 1.0 / E, alpha = 1.0, base = 10.0,
                          beta0 = -b1 * (0.75 / 0.09), beta1 = b1, V20 = V20, VT = V2.VT, clock_domain = V2.clock_domain,
                          appearance_domain = V2.appearance_domain, sink_domain = V2.sink_domain, vegetation_domain = V2.vegetation_domain))
    gpo = update_params(gpo, ["k_m_leaf", "k_m_stem", "k_m_fruit", "gen_leaf_target", "gen_N0", "set_r_low", "set_r_high",
                              "DMC_fruit", "harvest_W"], [0.03, 0.015, 0.01, 25.0, 4.0, 0.3, 1.0, 0.035, 323.0 * 0.035])
    ro = simulate_growth_measured(update_params(pco, ["J_max"], [2.0e-4]), gpo, seaso; LAI0 = 0.25, fruit_keep = keepo,
                                  generic_params = genp(60.7, 275.0, 0.0882, 3.2))
    # GreenhouseTwin
    CRJ = joinpath(EQ_HERE, "..", "configfiles", "crop_n105.json")
    gpn = GT.load_params(CRJ); pcn = pN
    O = GT.crop_observations(EQ_INV)
    zc = vcat([getproperty(gpn, Symbol(q)) for q in GT.CROP_QOI], gpn.sigma_y, gpn.sigma_leaf, gpn.sigma_fd, gpn.sigma_n)
    rn = GT.crop_run(zc, pcn, gpn, O)
    eqcheck("  measured season: same days and samples (T, CO2 identical; PAR to rounding)",
            O.season.ndays == seaso.ndays && all(size(a) == size(b) && a[1:2, :] == b[1:2, :] && maximum(abs.(a[3, :] .- b[3, :]); init = 0.0) < 1e-10
                                                for (a, b) in zip(O.season.day_samples, seaso.day_samples)))
    eqcheck("  kept-fruit share identical", O.keep == keepo)
    for (nm, a, b, tol) in (("yield FW", rn.yield_FW, ro.yield_FW, 1e-9), ("LAI", rn.LAI, ro.LAI, 1e-9), ("node", rn.node, ro.node, 1e-9),
                            ("leaves", rn.leaves, ro.leaves, 0.0), ("harvested fruits", rn.harv_n, ro.harv_n, 1e-9),
                            ("Pg", rn.Pg, ro.Pg, 1e-9), ("W_stem", rn.W_stem, ro.W_stem, 1e-9))
        d = maximum(abs.(a .- b)); eqcheck(@sprintf("  %-17s max |diff| %.2e (final %.4g)", nm, d, a[end]), d <= tol * max(1.0, maximum(abs.(b))))
    end
    hn = [isnan(x) ? -1.0 : x for x in rn.harv_age]; ho = [isnan(x) ? -1.0 : x for x in ro.harv_age]
    eqcheck("  harvest ages identical", maximum(abs.(hn .- ho)) < 1e-9)
    Un = GT.crop_energy(zc, pcn, gpn, O)
    @printf("    crop energy at the prior centre (GreenhouseTwin): %.6f\n", Un)
end

# ---------------------------------------------------------------------------------------------------------------------
println("\n== 9. rhs! timing and allocations (same state and time) ==")
let u = copy(SN.u0), du8 = zeros(8), du10 = zeros(10), d = 50
    co = GreenhouseSim.SimContext(pO, weather_c, ctrl_x, t0 + 86400.0 * d, 86400.0 * d)
    cn = GT.SimContext(pN, SN.weather, SN.controls, SN.t0 + 86400.0 * d, 86400.0 * d)
    uo = old_u(u)
    GreenhouseSim.rhs!(du10, uo, co, 43200.0); GT.rhs!(du8, u, cn, 43200.0)
    ao = @allocated GreenhouseSim.rhs!(du10, uo, co, 43200.0); an = @allocated GT.rhs!(du8, u, cn, 43200.0)
    bo = @belapsed GreenhouseSim.rhs!($du10, $uo, $co, 43200.0); bn = @belapsed GT.rhs!($du8, $u, $cn, 43200.0)
    eqcheck(@sprintf("  allocations per rhs! call: V2.2 %d, new %d", ao, an), an == 0)
    eqcheck(@sprintf("  time per rhs! call: V2.2 %.0f ns, new %.0f ns", 1e9bo, 1e9bn), bn <= 1.2 * bo)
end

@printf("\n%d PASS, %d FAIL\n", EQ_NP[], EQ_NF[])
