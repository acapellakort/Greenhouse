# =============================================================================
# check_solar_column.jl -- Note 103: energy-conserving solar column (flag solar_cons), LAI-dependent canopy optics
# =============================================================================
# Run from the GreenhouseSim folder (seconds, no data files):
#     julia --project=. test/check_solar_column.jl
# Calls the ENGINE kernels (solar_column, climate_rhs, floor_full_fluxes, cover_T, phys_par) at one fixed state,
# n79 optics, I2 = 500 W m-2; prints PASS/FAIL per check and the LAI / screen tables of Note 103.
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim, Printf
const GS  = GreenhouseSim
const CF  = joinpath(@__DIR__, "..", "..", "..", "invernaderos2.0", "code", "configfiles")
p0 = load_params(joinpath(CF, "constants_climate.json"), joinpath(CF, "constants_cropphoto.json"))
base = update_params(p0, ["tau1", "solar_split", "floor_full", "use_cover_node", "phys_light", "f_mass", "cover_area", "mass_node"],
                     [0.65 / 0.9, 1.0, 1.0, 1.0, 1.0, 0.15, 1.6, 1.0])
pON  = update_params(base, ["solar_cons", "solar_mass_share"], [1.0, 1.0])
pOFF = update_params(base, ["solar_cons"], [0.0])
npass = Ref(0); nfail = Ref(0)
check(name, ok) = (ok ? (npass[] += 1) : (nfail[] += 1); @printf("%-84s %s\n", name, ok ? "PASS" : "FAIL"))
LAIs  = (0.0, 0.5, 2.3, 5.0, 10.0)
cases = [("retracted", 0.0, 0.0), ("ES 50 %", 0.5, 0.0), ("ES 100 %", 1.0, 0.0), ("BL 50 %", 0.0, 0.5),
         ("BL 100 %", 0.0, 1.0), ("ES+BL 100 %", 1.0, 1.0)]
I2 = 500.0

println("\n== 1. canopy layer optics (t = exp(-k LAI), r = rho_inf (1 - t), a = 1 - t - r) ==")
for (bn, k, ri) in (("PAR", pON.beta1, pON.rho1), ("NIR", pON.beta3, pON.rho_can_nir))
    tprev = 2.0
    for L in LAIs
        t, r = canopy_layer(k, ri, L); a = 1 - t - r
        @printf("  %s LAI %5.1f: t %.4f  r %.4f  a %.4f\n", bn, L, t, r, a)
        L == 0 && check("  $bn LAI 0: t = 1, r = a = 0", t == 1 && r == 0 && a == 0)
        check(@sprintf("  %s LAI %.1f: t, r, a in [0,1]; t decreases with LAI", bn, L), all(0 .<= (t, r, a) .<= 1) && t < tprev)
        tprev = t
    end
    t, r = canopy_layer(k, ri, 40.0)
    check("  $bn LAI 40: direct transmission -> 0, reflection -> rho_inf", t < 1e-4 && abs(r - ri) < 1e-4)
end

println("\n== 2. compact canopy-floor formulas = column with transparent roof and screens ==")
for L in LAIs
    t, r = canopy_layer(pON.beta3, pON.rho_can_nir, L); a = 1 - t - r; f = 1 - pON.alpha7; D = 1 - r * f
    c = column4(1.0, 0.0, 1.0, 0.0, 1.0, 0.0, t, r, f)
    check(@sprintf("  NIR LAI %.1f: Qcan = aI(1+ft/D), Qflr = (1-f)tI/D, Rup = I(r+ft^2/D)", L),
          isapprox(c.A4, a * (1 + f * t / D); atol = 1e-14) && isapprox(c.Af, (1 - f) * t / D; atol = 1e-14) &&
          isapprox(c.Rout, r + f * t^2 / D; atol = 1e-14))
end

println("\n== 3. energy balance per band and in W m-2 (all LAI x all screen states) ==")
trp, rrp, trn, rrn = roof_props(pON)
@printf("  roof: tau_par %.3f rho_par %.3f (abs %.3f = eta1 %.2f + cov_abs %.2f)\n", trp, rrp, 1 - trp - rrp, pON.eta1, pON.cov_abs)
for L in LAIs, (nm, ue, ub) in cases
    for (bn, tr, rr, tes, res, tbl, rbl, k, ri, rf) in
        (("PAR", trp, rrp, pON.scr_es_tau_par, pON.scr_es_rho_par, pON.scr_bl_tau_par, pON.scr_bl_rho_par, pON.beta1, pON.rho1, pON.rho2),
         ("NIR", trn, rrn, pON.scr_es_tau_nir, pON.scr_es_rho_nir, pON.scr_bl_tau_nir, pON.scr_bl_rho_nir, pON.beta3, pON.rho_can_nir, 1 - pON.alpha7))
        es = layer_eff(tes, res, ue); bl = layer_eff(tbl, rbl, ub); ct = canopy_layer(k, ri, L)
        c = column4(tr, rr, es..., bl..., ct..., rf)
        v = (c.Rout, c.A1, c.A2, c.A3, c.A4, c.Af)
        (abs(sum(v) - 1) < 1e-12 && all(0 .<= v .<= 1)) ||
            check(@sprintf("  %s LAI %.1f %s: sum %.15f", bn, L, nm, sum(v)), false)
    end
    so = solar_column(I2, ue, ub, L, pON)
    s = so.Qcan + so.Qflr + so.Qair + so.Qcov + so.Qout
    abs(s - I2) < 1e-9 || check(@sprintf("  W balance LAI %.1f %s: %.9f", L, nm, s), false)
end
check("  per band: sky + roof + 2 screens + canopy + floor = 1, all in [0,1] (30 states x 2 bands)", true)
check("  W m-2: canopy + floor-level + air + cover + sky = I2 (30 states)", true)
let so = solar_column(0.0, 1.0, 1.0, 2.3, pON)
    check("  I2 = 0: every solar flux is 0", all(getfield(so, k) == 0 for k in fieldnames(SolarOut) if !(k in (:mode, :Spar, :Snir))))
end

println("\n== 4. engine books each flux once (cover node off isolates the solar terms) ==")
T1, T2, C1, V1, T5 = 295.0, 296.0, 700.0, 2000.0, 293.0
I3, I4, I5, I6, I8, I10, I11, I9p = 300.0, 270.0, 288.0, 273.15, 2.0, 834.7, 1200.0, 150.0
rhs(I2, L, U1, so, p) = climate_rhs(T1, T2, C1, V1, L, I2, I3, I4, I5, I6, T5, I8, I9p, I10, I11, 0.0,
                                    U1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, p; so = so)
ffl(I2, L, U1, so, p) = GS.floor_full_fluxes(T1, T2, T5, L, I2, I4, I5, U1, p; V1 = V1, so = so)
pn = update_params(pON, ["use_cover_node"], [0.0])
for L in (0.5, 2.3, 5.0), (nm, ue, ub) in cases
    so = solar_column(I2, ue, ub, L, pn); s0 = solar_column(0.0, ue, ub, L, pn)
    dcan = (rhs(I2, L, ue, so, pn)[1] - rhs(0.0, L, ue, s0, pn)[1]) * pn.alpha1 * L
    dair = (rhs(I2, L, ue, so, pn)[2] - rhs(0.0, L, ue, s0, pn)[2]) * pn.phi2 * pn.rho3 * pn.alpha5
    dflr = ffl(I2, L, ue, so, pn)[3] - ffl(0.0, L, ue, s0, pn)[3]
    ok = isapprox(dcan, so.Qcan; rtol = 1e-9) && isapprox(dair, so.Qair; rtol = 1e-9) &&
         isapprox(dflr, (1 - pn.f_mass) * so.Qflr; rtol = 1e-9)
    check(@sprintf("  LAI %.1f %-11s canopy = Qcan, air = Qair (no r8), floor = (1-f_mass) Qflr", L, nm), ok)
end
let so = solar_column(I2, 1.0, 0.0, 2.3, pON)
    kR = 4.0 * pON.epsil3 * pON.sigma * (0.5 * (T2 + I4))^3
    dTc = GS.cover_T(T2, I4, I5, I2, pON; Rx = so.Qcov, cons = true) - GS.cover_T(T2, I4, I5, 0.0, pON)
    check("  cover node gets Qcov once and no cov_abs I2 on top", isapprox(dTc, so.Qcov / (2kR + pON.h_cover_out); rtol = 1e-12))
    pp = GS.phys_par(999.0, 1.0, 0.0, (BlackScr = t -> 0.0,), 0.0, pON; so = so)
    check("  plant sun PAR = column ParTop (same field as the heat terms)", pp == so.ParTop)
end
check("  flag off: climate_rhs(so = nothing) is the V2.2 call", rhs(I2, 2.3, 0.5, nothing, pOFF) == rhs(I2, 2.3, 0.5, nothing, base))

let pb = update_params(p0, ["solar_cons"], [0.0])
    check("  rhs! path: SolarOut is concrete (isconcretetype)", isconcretetype(SolarOut) && solar_mode(SOLAR_NONE) == 0)
end

println("\n== 5. tables (I2 = 500 W m-2) ==")
println("  screens retracted, by LAI: canopy PAR | canopy NIR | floor-level PAR | floor-level NIR | sky | plant PAR on top")
for L in (0.0, 0.5, 2.3, 5.0)
    so = solar_column(I2, 0.0, 0.0, L, pON)
    @printf("  LAI %4.1f  %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n", L, so.QcanP, so.QcanN, so.QflrP, so.QflrN, so.Qout, so.ParTop)
end
println("  LAI 2.3, by screen state: plant PAR | canopy | floor-level | air | cover | screens | sky | sum")
for (nm, ue, ub) in cases
    so = solar_column(I2, ue, ub, 2.3, pON)
    @printf("  %-12s %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %8.3f\n", nm, so.ParTop, so.Qcan, so.Qflr, so.Qair, so.Qcov,
            so.Qscr, so.Qout, so.Qcan + so.Qflr + so.Qair + so.Qcov + so.Qout)
    nm == "retracted" || check(@sprintf("  %-11s deployment lowers plant PAR, canopy and floor together", nm),
        so.ParTop < solar_column(I2, 0.0, 0.0, 2.3, pON).ParTop && so.Qcan < solar_column(I2, 0.0, 0.0, 2.3, pON).Qcan &&
        so.Qflr < solar_column(I2, 0.0, 0.0, 2.3, pON).Qflr)
end
@printf("\n%d PASS, %d FAIL\n", npass[], nfail[])
