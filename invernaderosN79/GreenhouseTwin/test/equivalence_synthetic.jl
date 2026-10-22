# =============================================================================
# equivalence_synthetic.jl -- V2.2 (n79 flags) vs GreenhouseTwin on SYNTHETIC inputs (no data files; ~1 min)
# =============================================================================
#     julia --project=. test/equivalence_synthetic.jl          (INV=<invernaderos root> to override the default)
# Same inputs given to both codes: 4 days of smooth weather with screens, vents, lamps, CO2 valve, a pipe that is
# logged warm part of the time, LAI, crop respiration and fresh mass. Compares rhs! at random states and times and
# 4-day trajectories. Parameters: V2.2 JSONs + the n79 fixed list + prior centres (as run_climate_calibration_cl.jl).
# =============================================================================
const SY_INV = get(ENV, "INV", normpath(joinpath(@__DIR__, "..", "..", "..")))
const SY_OLD = joinpath(SY_INV, "invernaderosV2.2", "GreenhouseSim", "src", "GreenhouseSim.jl")
include(SY_OLD); include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
const GS = GreenhouseSim; const GT = GreenhouseTwin
using Printf, Random, OrdinaryDiffEq

# --- parameters -------------------------------------------------------------------------------------------------------
CF = joinpath(SY_INV, "invernaderos2.0", "code", "configfiles")
pO = GS.load_params(joinpath(CF, "constants_climate.json"), joinpath(CF, "constants_cropphoto.json"))
FIX = ["nu4" => 4.46e-5, "U8_min_leak" => 0.0, "heat_ctrl" => 1.0, "hc_Kp" => 10.0, "hc_Ki" => 10.0 / 1800, "vent_screen" => 1.0,
       "stomata_fix" => 1.0, "vapour_fix" => 1.0, "cond_heat" => 1.0, "f_cond" => 0.0, "vent_ctrl" => 0.0, "vc_hband" => 1.0,
       "pipe_lag" => 0.0, "black_scr" => 0.0, "hc_log_pipe" => 1.0, "vc_min_log" => 1.0, "soil_dyn" => 1.0, "floor_full" => 1.0,
       "solar_split" => 1.0, "deep_k" => 1.0, "mass_node" => 1.0, "sd_tau" => 14 * 86400.0, "lamp_full" => 1.0, "vent_split" => 1.0,
       "vent_exp" => 1.0, "vs_clee" => 1.0, "phi2" => 5.4, "alpha5" => 1005.0, "nu5" => 3000.0, "gamma1" => 1.0, "n_pipes" => 1.0,
       "tau1" => 0.65 / 0.9, "psi2" => 35556.0, "cover_area" => 1.6, "cover_cond" => 1.0, "co2_resp" => 1.0, "hc_dT" => 0.0,
       "lamp_scr_loss" => 0.0, "s_rad" => 0.0, "dew_bl" => 1.0, "leak_shared" => 1.0, "rd_q10" => 1.0, "leak_stack" => 1.0,
       "lamp_Pel" => 117.0, "crop_cap" => 1.0, "crop_cp" => 3800.0, "solar_cons" => 1.0, "solar_mass_share" => 1.0, "par_below" => 1.0]
pO = GS.update_params(pO, first.(FIX), last.(FIX))
pN = GT.load_params(joinpath(@__DIR__, "..", "configfiles", "climate_n79.json"))
pO = GS.update_params(pO, GT.CLIMATE_QOI, [getproperty(pN, Symbol(q)) for q in GT.CLIMATE_QOI])

# --- synthetic inputs (5-min, 4 days) ------------------------------------------------------------------------------------
const T0 = 1.5342e9
n = 4 * 288 + 1; ts = [T0 + 300.0 * (k - 1) for k in 1:n]; tr = ts .- T0; h = (tr ./ 3600) .% 24
sun = max.(0.0, sin.(pi .* (h .- 6) ./ 13))
Ig = 650 .* sun .* (0.8 .+ 0.2 .* cos.(tr ./ 5000)); Tout = 285 .+ 7 .* sun .+ 0.5 .* sin.(tr ./ 9000)
Tsky = Tout .- 12 .+ 8 .* (sin.(tr ./ 40000) .> 0); Wind = 2.5 .+ 1.5 .* sin.(tr ./ 7000) .^ 2
VPo = 0.75 .* GS.Pws.(Tout); LAI = 1.5 .+ tr ./ (4 * 86400)
U1 = clamp.(1 .- 3 .* sun, 0.0, 1.0); Bl = Float64.((h .> 20) .| (h .< 1)) .* 0.8
Lee = clamp.(2 .* sun .- 0.5, 0.0, 1.0); Wv = 0.5 .* Lee; U10 = Float64.(sun .> 0.2) .* 0.6
U11 = fill(291.15, n); U12 = Float64.((h .> 2) .& (h .< 8)); Tpipe = [isodd(floor(Int, x / 4)) ? 318.0 : Tout[k] + 2.0 for (k, x) in enumerate(h)]
CR = fill(0.004, n); FW = 3.0 .+ tr ./ 86400 .* 0.05
A(v) = (let tt = ts, vv = v; t -> vv[clamp(searchsortedlast(tt, t), 1, length(vv))] end)
Z(v) = (let tt = tr, vv = v; t -> vv[clamp(searchsortedlast(tt, t), 1, length(vv))] end)
zc = Z(zeros(n))
wO = GS.WeatherInputs(A(Ig), A(Tsky), A(Tout), _t -> 273.15, _t -> 291.15, A(Wind), _t -> 834.7, A(VPo), A(Ig), A(LAI))
cO = GS.ControlInputsX(Z(U1), zc, zc, zc, zc, zc, Z(Wv), Z(Lee), zc, Z(U10), Z(U11), Z(U12), Z(Tpipe),
                       Z(fill(298.0, n)), Z(fill(3.0, n)), Z(Bl), Z(zeros(n)), Z(CR), Z(FW))
SF(t, v) = GT.StepFun(t, v); zN = SF(tr, zeros(n))
wN = GT.WeatherInputs(SF(ts, Ig), SF(ts, Tsky), SF(ts, Tout), SF(ts, Wind), SF(ts, VPo), GT.ConstFun(834.7), GT.ConstFun(273.15), SF(ts, LAI))
cN = GT.ControlInputs(SF(tr, U1), SF(tr, Bl), zN, zN, zN, zN, zN, zN, SF(tr, Lee), SF(tr, Wv), zN, SF(tr, U10), SF(tr, U11),
                      SF(tr, U12), SF(tr, Tpipe), SF(tr, CR), SF(tr, FW))
const IMAP = [1, 2, 3, 4, 5, 6, 9, 10]
old_u(u) = [u[1], u[2], u[3], u[4], u[5], u[6], 0.0, u[2], u[7], u[8]]

np = Ref(0); nf = Ref(0)
chk(name, ok) = (ok ? (np[] += 1) : (nf[] += 1); @printf("%-88s %s\n", name, ok ? "PASS" : "FAIL"))

println("== rhs! at 5000 random states and times ==")
let rng = MersenneTwister(1), N = 5000, DO = zeros(N, 8), DN = zeros(N, 8), duo = zeros(10), dun = zeros(8)
    for k in 1:N
        d = rand(rng, 0:2); t = isodd(k) ? 300.0 * rand(rng, 0:287) : 86400.0 * rand(rng)
        T2 = 285 + 20rand(rng); RH = 40 + 58rand(rng)
        u = [T2 + 3randn(rng), T2, 400 + 1100rand(rng), RH / 100 * GS.Pws(T2), 3000rand(rng) * (rand(rng) < 0.3 ? -1 : 1),
             T2 - 5 + 8rand(rng), 285 + 10rand(rng), T2 - 3 + 6rand(rng)]
        GS.rhs!(duo, old_u(u), GS.SimContext(pO, wO, cO, T0 + 86400.0 * d, 86400.0 * d), t)
        GT.rhs!(dun, u, GT.SimContext(pN, wN, cN, T0 + 86400.0 * d, 86400.0 * d), t)
        DO[k, :] = duo[IMAP]; DN[k, :] = dun
    end
    for (j, nm) in enumerate(("dT1", "dT2", "dC1", "dV1", "de_int", "dT5", "dTd", "dTm"))
        sc = maximum(abs.(DO[:, j])); e = maximum(abs.(DO[:, j] .- DN[:, j])) / sc
        chk(@sprintf("  %-7s max |diff| / max |value| %.2e (scale %.3g)", nm, e, sc), e < 1e-10)
    end
end

println("== 3-day trajectories (judged against the solver's own sensitivity) ==")
# The trajectories react to rounding-level differences (discontinuous inputs, PI and max() switches, adaptive steps):
# V2.2 against itself with the initial state scaled by (1 + 1e-14) diverges by the same order. A trajectory check
# passes when |new - V2.2| is within 5x that noise floor (plus 1e-9 of the range).
function run3(f, ctxf, u0, rt)
    u = copy(u0); out = Matrix{Float64}(undef, length(u0), 0)
    for d in 0:2
        s = solve(ODEProblem(f, u, (0.0, 86400.0), ctxf(d)), Tsit5(); saveat = 300.0, reltol = rt)
        out = hcat(out, Array(s)); u = s[:, end]
    end
    out
end
ctxO(d) = GS.SimContext(pO, wO, cO, T0 + 86400.0 * d, 86400.0 * d)
ctxN(d) = GT.SimContext(pN, wN, cN, T0 + 86400.0 * d, 86400.0 * d)
let u0 = [292.0, 292.5, 700.0, 1600.0, 0.0, 291.0, 290.0, 292.0]
    for rt in (1e-6, 1e-9)
        uo = run3(GS.rhs!, ctxO, old_u(u0), rt); up = run3(GS.rhs!, ctxO, old_u(u0 .* (1 + 1e-14)), rt); un = run3(GT.rhs!, ctxN, u0, rt)
        for (j, nm) in ((1, "T1 [K]"), (2, "T2 [K]"), (3, "C1 [mg m-3]"), (4, "V1 [Pa]"), (6, "T5 [K]"), (7, "Td [K]"), (8, "Tm [K]"))
            e = maximum(abs.(uo[IMAP[j], :] .- un[j, :])); b = maximum(abs.(uo[IMAP[j], :] .- up[IMAP[j], :]))
            rng_ = -(reverse(extrema(un[j, :]))...)
            chk(@sprintf("  reltol %.0e %-12s |new - V2.2| %.2e   noise floor |V2.2 - V2.2(u0+1e-14)| %.2e", rt, nm, e, b),
                e <= 5b + 1e-9 * abs(rng_))
        end
    end
end

println("== solar ==")
let so = GS.solar_column(500.0, 0.7, 0.3, 2.3, pO), sn = GT.solar_column(500.0, 0.7, 0.3, 2.3, pN)
    chk("  solar fluxes identical (all fields)", all(getfield(sn, f) == getfield(so, f) for f in fieldnames(GT.SolarFluxes)))
    chk("  both-faces factor identical", GS.par_both_faces_factor(2.3, pO) == GT.par_both_faces_factor(2.3, pN))
end

println("== crop: 116 synthetic days (T, CO2, PAR), n105 set-up at the prior centre ==")
let rng = MersenneTwister(3), nd = 116
    samples = [let m = 288; M = zeros(3, m)
                   for j in 1:m
                       hh = (j - 1) / 12; sl = max(0.0, sin(pi * (hh - 6) / 13))
                       M[1, j] = 291 + 6sl + randn(rng); M[2, j] = 700 + 300rand(rng); M[3, j] = 120sl * (0.5 + rand(rng)) + 40.65 * (hh > 2 && hh < 8)
                   end; M end for _ in 1:nd]
    keep = [0.5 + 0.17 * isodd(d ÷ 7) for d in 0:nd-1]
    # V2.2
    pco = GS.update_params(GS.load_params(joinpath(CF, "constants_climate.json"), joinpath(CF, "constants_cropphoto.json")),
                           ["par_below", "tau1", "J_max"], [1.0, 0.65 / 0.9, 2.0e-4])
    gpo = GS.load_growth_params(joinpath(SY_INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json"))
    gpo = GS.update_params(gpo, ["set_start_day", "harvest_W", "gen_clampT", "gen_parlight", "gen_beta", "gen_kappa", "gen_A",
                                 "gen_V", "gen_zclamp", "stem_density", "gen_canopy", "gen_doy0", "gen_fL", "gen_flower_node",
                                 "k_m_leaf", "k_m_stem", "k_m_fruit", "gen_leaf_target", "gen_N0", "set_r_low", "set_r_high", "DMC_fruit"],
                           [15.0, 323.0 * 0.035, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.6, 1.0, 226.0, 0.59, 5.0,
                            0.03, 0.015, 0.01, 25.0, 4.0, 0.3, 1.0, 0.035])
    V2 = GS.CUCUMBER_PARAMS_V2; AT = 60.7 * V2.s * (V2.p - 1) / 275.0
    G2 = GS.GenericParameters(s = V2.s, p = V2.p, peak = V2.peak, A0 = -10.0 * AT, AT = AT, kappa = 1.0 / 275.0, alpha = 1.0,
                              base = 10.0, beta0 = -0.0882 * (0.75 / 0.09), beta1 = 0.0882, V20 = 3.2, VT = V2.VT,
                              clock_domain = V2.clock_domain, appearance_domain = V2.appearance_domain,
                              sink_domain = V2.sink_domain, vegetation_domain = V2.vegetation_domain)
    ro = GS.simulate_growth_measured(pco, gpo, GS.MeasuredSeason(nd, samples, 300.0); LAI0 = 0.25, fruit_keep = keep, generic_params = G2)
    # GreenhouseTwin
    gpn = GT.load_params(joinpath(@__DIR__, "..", "configfiles", "crop_n105.json"))
    Gn = GT.cucumber_generic(gpn.W_pot, gpn.E_fruit, gpn.app_beta1, gpn.V20, gpn)
    chk("  generic parameters identical (s, p, peak, A0, AT, kappa, beta0, beta1, V20, VT)",
        (Gn.s, Gn.p, Gn.peak, Gn.A0, Gn.AT, Gn.kappa, Gn.beta0, Gn.beta1, Gn.V20, Gn.VT) ==
        (G2.s, G2.p, G2.peak, G2.A0, G2.AT, G2.kappa, G2.beta0, G2.beta1, G2.V20, G2.VT))
    rn = GT.simulate_crop_measured(GT.update_params(pN, ["J_max"], [2.0e-4]), gpn, Gn, GT.MeasuredSeason(nd, samples, 300.0), keep)
    for (nm, a, b) in (("yield FW", rn.yield_FW, ro.yield_FW), ("LAI", rn.LAI, ro.LAI), ("node", rn.node, ro.node),
                       ("leaves", rn.leaves, ro.leaves), ("harvested n", rn.harv_n, ro.harv_n), ("Pg", rn.Pg, ro.Pg),
                       ("W_leaf", rn.W_leaf, ro.W_leaf), ("W_stem", rn.W_stem, ro.W_stem), ("W_fruit", rn.W_fruit, ro.W_fruit))
        e = maximum(abs.(a .- b)); chk(@sprintf("  %-12s max |diff| %.2e (end %.4g)", nm, e, a[end]), e <= 1e-12 * max(1.0, maximum(abs.(b))))
    end
    ha(x) = [isnan(v) ? -1.0 : v for v in x]
    chk("  harvest ages identical", ha(rn.harv_age) == ha(ro.harv_age))
end
@printf("\n%d PASS, %d FAIL\n", np[], nf[])
