# =============================================================================
# check_solar_column.jl -- solar column: optics, energy closure, booking in the model (seconds, no data files)
# =============================================================================
#     julia --project=. test/check_solar_column.jl
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin, Printf
const GT = GreenhouseTwin
p = load_params(joinpath(@__DIR__, "..", "configfiles", "climate_n79.json"))
npass = Ref(0); nfail = Ref(0)
check(name, ok) = (ok ? (npass[] += 1) : (nfail[] += 1); @printf("%-90s %s\n", name, ok ? "PASS" : "FAIL"))
LAIs  = (0.0, 0.5, 2.3, 5.0, 10.0)
cases = [("retracted", 0.0, 0.0), ("ES 50 %", 0.5, 0.0), ("ES 100 %", 1.0, 0.0), ("BL 50 %", 0.0, 0.5),
         ("BL 100 %", 0.0, 1.0), ("ES+BL 100 %", 1.0, 1.0)]
I2 = 500.0

println("\n== 1. canopy layer (t = exp(-k LAI), r = rho_inf (1 - t)) ==")
for (bn, k, ri) in (("PAR", p.beta1, p.rho1), ("NIR", p.beta3, p.rho_can_nir))
    tprev = 2.0
    for L in LAIs
        t, r = canopy_layer(k, ri, L); a = 1 - t - r
        L == 0 && check("  $bn LAI 0: t = 1, r = a = 0", t == 1 && r == 0 && a == 0)
        check(@sprintf("  %s LAI %.1f: t, r, a in [0,1]; t decreases with LAI", bn, L), all(0 .<= (t, r, a) .<= 1) && t < tprev)
        tprev = t
    end
    t, r = canopy_layer(k, ri, 40.0)
    check("  $bn LAI 40: transmission -> 0, reflection -> rho_inf", t < 1e-4 && abs(r - ri) < 1e-4)
end

println("\n== 2. closure per band and in W m-2 (all LAI x screen states) ==")
trp, rrp, _, _ = roof_props(p)
@printf("  roof: tau %.3f rho %.3f absorption %.3f = eta1 %.2f + cov_abs %.2f\n", trp, rrp, 1 - trp - rrp, p.eta1, p.cov_abs)
check("  roof absorption = eta1 + cov_abs", isapprox(1 - trp - rrp, p.eta1 + p.cov_abs; atol = 1e-14))
okb = true; okw = true
for L in LAIs, (nm, ue, ub) in cases
    for (tes, res, tbl, rbl, k, ri, rf) in ((p.scr_es_tau_par, p.scr_es_rho_par, p.scr_bl_tau_par, p.scr_bl_rho_par, p.beta1, p.rho1, p.rho2),
                                             (p.scr_es_tau_nir, p.scr_es_rho_nir, p.scr_bl_tau_nir, p.scr_bl_rho_nir, p.beta3, p.rho_can_nir, 1 - p.alpha7))
        c = column4(trp, rrp, layer_eff(tes, res, ue)..., layer_eff(tbl, rbl, ub)..., canopy_layer(k, ri, L)..., rf)
        v = (c.Rout, c.A1, c.A2, c.A3, c.A4, c.Af)
        global okb &= abs(sum(v) - 1) < 1e-12 && all(0 .<= v .<= 1)
    end
    s = solar_column(I2, ue, ub, L, p)
    global okw &= abs(s.Qcan + s.Qflr + s.Qair + s.Qcov + s.Qout - I2) < 1e-9
end
check("  per band: sky + roof + 2 screens + canopy + floor = 1 (30 states x 2 bands)", okb)
check("  W m-2: canopy + floor level + air + cover + sky = Iglob (30 states)", okw)

println("\n== 3. both-faces PAR and the crop-path factor ==")
for L in (0.5, 2.3, 5.0), (nm, ue, ub) in cases
    s = solar_column(I2, ue, ub, L, p); t, r = canopy_layer(p.beta1, p.rho1, L); a = 1 - t - r
    check(@sprintf("  LAI %.1f %-11s ParTop + ParBot = absorbed canopy PAR / a", L, nm),
          isapprox(s.ParTop + s.ParBot, s.QcanP / a; rtol = 1e-12) && s.ParBot >= 0)
end
for L in (0.0, 0.5, 2.3, 5.0)
    s0 = solar_column(1.0, 0.0, 0.0, L, p)
    g = par_both_faces_factor(L, p)
    check(@sprintf("  g(LAI %.1f) = %.3f = (ParTop + ParBot)/(eta2 tau_roof) >= 1", L, g),
          isapprox(g, (s0.ParTop + s0.ParBot) / (p.eta2 * trp); rtol = 1e-12) && g >= 1)
end

println("\n== 4. tables (Iglob 500 W m-2) ==")
println("  LAI 2.3: plant PAR top | canopy | floor level | air | cover | screens | sky | sum")
for (nm, ue, ub) in cases
    s = solar_column(I2, ue, ub, 2.3, p)
    @printf("  %-12s %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %8.3f\n", nm, s.ParTop, s.Qcan, s.Qflr, s.Qair, s.Qcov, s.Qscr,
            s.Qout, s.Qcan + s.Qflr + s.Qair + s.Qcov + s.Qout)
end
@printf("\n%d PASS, %d FAIL\n", npass[], nfail[])
