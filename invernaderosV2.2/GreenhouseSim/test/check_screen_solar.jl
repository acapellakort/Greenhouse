# =============================================================================
# check_screen_solar.jl -- Note 101: screen optics shared by plant light and solar heat (flag scr_solar)
# =============================================================================
# Run from the GreenhouseSim folder (seconds, no data files needed):
#     julia --project=. test/check_screen_solar.jl
# Calls the ENGINE kernels (screen_solar, climate_rhs, floor_full_fluxes, phys_par) at one fixed state, n79 optics
# (tau1 = 0.65/0.9, solar_split, floor_full, cover node), I2 = 500 W m-2, LAI 2.2, and prints PASS/FAIL per check
# plus the before/after table of Note 101.
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim, Printf
const GS  = GreenhouseSim
const INV = joinpath(@__DIR__, "..", "..", "..")
const CF  = joinpath(INV, "invernaderos2.0", "code", "configfiles")
p0 = load_params(joinpath(CF, "constants_climate.json"), joinpath(CF, "constants_cropphoto.json"))
base = update_params(p0, ["tau1", "solar_split", "floor_full", "use_cover_node", "phys_light", "f_mass", "cover_area"],
                     [0.65 / 0.9, 1.0, 1.0, 1.0, 1.0, 0.15, 1.6])
pOFF = update_params(base, ["scr_solar"], [0.0]); pON = update_params(base, ["scr_solar"], [1.0])

# fixed state and inputs (identical in every call)
T1, T2, C1, V1, T5 = 295.0, 296.0, 700.0, 2000.0, 293.0
LAI, I3, I4, I5, I6, I8, I10, I11 = 2.2, 300.0, 270.0, 288.0, 273.15, 2.0, 834.7, 1200.0
I9phys = 150.0                       # physiology PAR passed to the stomata: SAME in on/off calls (isolates heat terms)
U12 = 0.0; A = 0.0
rhs(I2, U1, so, p) = climate_rhs(T1, T2, C1, V1, LAI, I2, I3, I4, I5, I6, T5, I8, I9phys, I10, I11, A,
                                 U1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, U12, p; so = so)
ffl(I2, U1, so, p) = GS.floor_full_fluxes(T1, T2, T5, LAI, I2, I4, I5, U1, p; V1 = V1, so = so)
heatC(p) = p.alpha1 * LAI                                  # J m-2 K-1 canopy
heatA(p) = p.phi2 * p.rho3 * p.alpha5                      # J m-2 K-1 air
npass = Ref(0); nfail = Ref(0)
check(name, ok) = (ok ? (npass[] += 1) : (nfail[] += 1); @printf("%-72s %s\n", name, ok ? "PASS" : "FAIL"))

cases = [("retracted", 0.0, 0.0), ("ES 50 %", 0.5, 0.0), ("ES 100 %", 1.0, 0.0), ("BL 50 %", 0.0, 0.5),
         ("BL 100 %", 0.0, 1.0), ("ES+BL 50 %", 0.5, 0.5), ("ES+BL 100 %", 1.0, 1.0)]

println("\n== 1. optics: bounds and conservation per band ==")
for (nm, ue, ub) in cases, (bnm, B) in (("PAR", :P), ("NIR", :N))
    f = getfield(screen_factors(ue, ub, pON), B)
    rp = roof_props(pON)                       # Note 103: roof rho derived (< 0 in the JSON)
    trb = B == :P ? rp[1] : rp[3]; rrb = B == :P ? rp[2] : rp[4]
    T, R = f.T, f.dR + rrb; Ar = f.dAr + (1 - trb - rrb)
    s = T + R + Ar + f.Aes + f.Abl
    check(@sprintf("%-12s %s  T+R+A_roof+A_es+A_bl = %.15f, all in [0,1]", nm, bnm, s),
          abs(s - 1) < 1e-12 && all(0 .<= (T, R, Ar, f.Aes, f.Abl) .<= 1))
end
let ok = true                        # random layers: conservation of the adding method
    for _ in 1:20000
        tr = rand(); rr = rand() * (1 - tr); t1 = rand(); r1 = rand() * (1 - t1); t2 = rand(); r2 = rand() * (1 - t2)
        a = layer_eff(t1, r1, rand()); b = layer_eff(t2, r2, rand()); o = stack3(tr, rr, a..., b...)
        ok &= abs(sum(o) - 1) < 1e-12 && all(0 .<= o .<= 1)
    end
    check("random stacks (20 000): sum = 1, bounded", ok)
end
let tr = 0.65, rr = 0.13, t = 0.75, r = 0.15
    check("two layers = Vanthoor tau12 = t1 t2 / (1 - r1 r2)", abs(stack3(tr, rr, t, r, 1.0, 0.0)[1] - tr * t / (1 - rr * r)) < 1e-14)
end

println("\n== 2. retracted screens: no extra attenuation (engine on == off, bitwise) ==")
so0 = screen_solar(500.0, 0.0, 0.0, pON)
check("S_par = S_nir = 1 and Qscr = Qcov = Qair = 0", so0.Spar == 1 && so0.Snir == 1 && so0.Qscr == 0 && so0.Qcov == 0)
check("climate_rhs(scr_solar on, retracted) == flag off", rhs(500.0, 0.0, so0, pON) == rhs(500.0, 0.0, nothing, pOFF))
check("floor_full_fluxes(on, retracted) == off", ffl(500.0, 0.0, so0, pON) == ffl(500.0, 0.0, nothing, pOFF))
c_ret = (BlackScr = t -> 0.0,)
check("phys_par(on, retracted) == I9 (sun)", GS.phys_par(162.5, 0.0, 0.0, c_ret, 0.0, pON; so = so0) == 162.5)

println("\n== 3. zero outside radiation: zero solar-derived gains (all screen states) ==")
for (nm, ue, ub) in cases
    so = screen_solar(0.0, ue, ub, pON)
    check(@sprintf("%-12s I2 = 0: Qscr = Qcov = Qair = Qrefl = 0; rhs on == off", nm),
          so.Qscr == 0 && so.Qcov == 0 && so.Qair == 0 && so.Qrefl == 0 &&
          rhs(0.0, ue, so, pON) == rhs(0.0, ue, nothing, pOFF) && ffl(0.0, ue, so, pON) == ffl(0.0, ue, nothing, pOFF))
end

println("\n== 4. deployment changes plant PAR and absorbed solar together; heat terms counted once ==")
p = pON; L = LAI; I2 = 500.0
I9sun = (1 - p.eta1) * p.tau1 * p.eta2 * I2
@printf("%-12s %6s %6s | %8s %8s | %8s %8s | %7s %7s | %7s %7s | %6s %6s %6s %6s\n", "case", "S_par", "S_nir",
        "PARold", "PARnew", "canOLD", "canNEW", "flrOLD", "flrNEW", "airOLD", "airNEW", "scrAbs", "toCov", "toAir", "refl")
for (nm, ue, ub) in cases
    so  = screen_solar(I2, ue, ub, p)
    # canopy solar from the engine: d(dT1) between I2 = 500 and I2 = 0 at the same screens, cover node OFF so the
    # cover temperature (and the long-wave terms) cannot change -> pure solar terms r2 + r3 + r5
    pn = update_params(p, ["use_cover_node"], [0.0]); pf = update_params(pOFF, ["use_cover_node"], [0.0])
    soN = screen_solar(I2, ue, ub, pn); soZ = screen_solar(0.0, ue, ub, pn)
    canN = (rhs(I2, ue, soN, pn)[1] - rhs(0.0, ue, soZ, pn)[1]) * heatC(pn)
    canO = (rhs(I2, ue, nothing, pf)[1] - rhs(0.0, ue, nothing, pf)[1]) * heatC(pf)
    airN = (rhs(I2, ue, soN, pn)[2] - rhs(0.0, ue, soZ, pn)[2]) * heatA(pn)       # r8 + screen share to the air
    airO = (rhs(I2, ue, nothing, pf)[2] - rhs(0.0, ue, nothing, pf)[2]) * heatA(pf)
    flN = ffl(I2, ue, soN, pn)[3] - ffl(0.0, ue, soZ, pn)[3]; flO = ffl(I2, ue, nothing, pf)[3] - ffl(0.0, ue, nothing, pf)[3]
    parO = I9sun * (1 - 0.98 * ub) * (1 - 0.25 * ue); parN = GS.phys_par(I9sun, ue, 0.0, (BlackScr = t -> ub,), 0.0, p; so = so)
    @printf("%-12s %6.3f %6.3f | %8.1f %8.1f | %8.1f %8.1f | %7.1f %7.1f | %7.1f %7.1f | %6.1f %6.1f %6.1f %6.1f\n",
            nm, so.Spar, so.Snir, parO, parN, canO, canN, flO, flN, airO, airN, so.Qscr, so.Qcov, so.Qair, so.Qrefl)
    # expected canopy change = analytic r2 + r3 (PAR, x S_par) + r5 sun (NIR, x S_nir): no other heat from PAR
    r4 = I9sun * so.Spar
    exp_can = r4 * (1 - p.rho1) * (1 - exp(-p.beta1 * L)) + r4 * exp(-p.beta1 * L) * p.rho2 * (1 - p.rho1) * (1 - exp(-p.beta2 * L)) +
              (1 - p.eta1) * p.alpha2 * p.eta3 * I2 * so.Snir
    exp_air = I2 * p.eta1 * (p.tau1 * p.eta2 * so.Spar + (p.alpha2 + p.alpha7) * p.eta3 * so.Snir) + so.Qair
    check(@sprintf("  %-10s canopy solar = r2 + r3 + r5(sun) x S only (PAR heat counted once)", nm), isapprox(canN, exp_can; rtol = 1e-9))
    check(@sprintf("  %-10s air solar = r8 x S + scr_abs_air x Qscr", nm), isapprox(airN, exp_air; rtol = 1e-9))
    check(@sprintf("  %-10s plant sun PAR = I9 x S_par (same factor as the heat terms)", nm), isapprox(parN, I9sun * so.Spar; rtol = 1e-12))
    ue + ub > 0 && check(@sprintf("  %-10s deployed: plant PAR and canopy/floor solar both reduced", nm),
                         parN < I9sun && canN < canO && flN < flO)
    # optics balance in W m-2: intercepted below-screen power = extra reflection + extra roof abs + screen abs
    Tb = p.eta2 * I2 * screen_factors(ue, ub, p).P.T + p.eta3 * I2 * screen_factors(ue, ub, p).N.T
    T0 = (1 - p.eta1) * p.tau1 * I2
    check(@sprintf("  %-10s W m-2: (T0 - T) = refl + roof + screens (%.2f)", nm, T0 - Tb),
          isapprox(T0 - Tb, so.Qrefl + so.Qroof + so.Qscr; atol = 1e-9))
end
let so = screen_solar(I2, 1.0, 1.0, pON)
    kR = 4.0 * pON.epsil3 * pON.sigma * (0.5 * (T2 + I4))^3
    dTc = GS.cover_T(T2, I4, I5, I2, pON; Rx = so.Qcov) - GS.cover_T(T2, I4, I5, I2, pON)
    check(@sprintf("cover node receives Qcov once: dTc = Qcov/(2kR + h_out) = %.3f K", dTc),
          isapprox(dTc, so.Qcov / (2kR + pON.h_cover_out); rtol = 1e-12))
end
println("\nPre-existing (not changed here): booked solar below the roof with screens retracted, LAI 2.2:")
let so = screen_solar(I2, 0.0, 0.0, pON), pf = update_params(pOFF, ["use_cover_node"], [0.0])
    can = (rhs(I2, 0.0, nothing, pf)[1] - rhs(0.0, 0.0, nothing, pf)[1]) * heatC(pf)
    air = (rhs(I2, 0.0, nothing, pf)[2] - rhs(0.0, 0.0, nothing, pf)[2]) * heatA(pf)
    flr = ffl(I2, 0.0, nothing, pf)[3] - ffl(0.0, 0.0, nothing, pf)[3]
    @printf("  transmitted (1-eta1) tau1 I2 = %.1f; canopy %.1f + floor %.1f + air r8 %.1f + mass f_mass I2 %.1f = %.1f W m-2\n",
            (1 - pOFF.eta1) * pOFF.tau1 * I2, can, flr, air, pOFF.f_mass * I2, can + flr + air + pOFF.f_mass * I2)
end
@printf("\n%d PASS, %d FAIL\n", npass[], nfail[])
