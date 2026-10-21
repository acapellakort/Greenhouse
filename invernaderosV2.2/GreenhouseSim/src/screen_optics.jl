# =============================================================================
# screen_optics.jl -- Note 101 (flag scr_solar, default off): one solar optical calculation for roof + screens,
# shared by plant light and solar heating.
# =============================================================================
# Bands: PAR and solar NIR (eta2, eta3 of the global radiation). Thermal LWIR (g2, tau3, r7, r11, r12, r15) is NOT
# handled here. Lamps hang below the screens and are not attenuated here.
#
# Layers, top to bottom: roof (the existing effective roof: PAR transmission (1 - eta1) tau1, i.e. 0.65 in n79),
# upper screen, lower screen. Each layer is symmetric (same tau, rho from both sides). A screen deployed by the
# fraction u in [0,1] (logged percentage / 100) has, as in Vanthoor (2011) eqs 8.14-8.17,
#     tau_eff = 1 - u (1 - tau),   rho_eff = u rho,   alpha_eff = u alpha = u (1 - tau - rho).
# Layers are combined with repeated reflections (adding method; for two layers tau12 = tau1 tau2 / (1 - rho1 rho2)).
# Approximation: radiation reflected upward by the canopy and the floor is not traced back to the screens.
#
# Output: screen factors S_b = T_b(screens) / T_b(retracted) by band. Every existing solar term below the screens
# is multiplied by its band's S_b, once, whatever roof convention the term carries (S_b = 1 exactly when the
# screens are retracted). The intercepted part is split into extra reflection (leaves through the roof), extra roof
# absorption, and screen absorption; screen absorption goes scr_abs_air to the crop air and the rest, with the extra
# roof absorption, to the algebraic cover/structure node.
# =============================================================================

"Screen deployed by fraction u: effective (tau, rho) (Vanthoor 2011, eqs 8.14-8.17)."
@inline layer_eff(tau, rho, u) = (1.0 - u * (1.0 - tau), u * rho)

"""
    stack3(tr, rr, t1, r1, t2, r2) -> (T, R, Ar, A1, A2)

Unit flux from above on three symmetric layers (roof tr, rr; upper screen t1, r1; lower screen t2, r2), all
repeated reflections included, nothing coming back from below. T + R + Ar + A1 + A2 = 1.
"""
@inline function stack3(tr, rr, t1, r1, t2, r2)
    d2  = 1.0 - r1 * r2
    r12 = r1 + t1 * t1 * r2 / d2            # screen pair seen from above
    D1  = tr / (1.0 - rr * r12)             # down, below the roof
    U1  = r12 * D1                          # up, below the roof
    D2  = t1 * D1 / d2                      # down, between the screens
    U2  = r2 * D2                           # up, between the screens
    T   = t2 * D2
    R   = rr + tr * U1
    Ar  = (1.0 - tr - rr) * (1.0 + U1)
    A1  = (1.0 - t1 - r1) * (D1 + U2)
    A2  = (1.0 - t2 - r2) * D2
    return T, R, Ar, A1, A2
end

"""
    roof_props(p) -> (tau_par, rho_par, tau_nir, rho_nir)

Roof (glass + structure) optics. tau_par = (1 - eta1) tau1 (0.65 in n79); tau_nir = roof_tau_nir, or tau_par when < 0.
Note 103: roof_rho_* < 0 -> derived so that the roof absorbs exactly what V2.2 already books there with the screens
retracted: alpha_roof = eta1 (building elements) + cov_abs (glass), rho = 1 - tau - eta1 - cov_abs (0.15 in n79).
"""
@inline function roof_props(p)
    trp = (1.0 - p.eta1) * p.tau1
    trn = p.roof_tau_nir < 0.0 ? trp : p.roof_tau_nir
    rrp = p.roof_rho_par < 0.0 ? max(1.0 - trp - p.eta1 - p.cov_abs, 0.0) : p.roof_rho_par
    rrn = p.roof_rho_nir < 0.0 ? max(1.0 - trn - p.eta1 - p.cov_abs, 0.0) : p.roof_rho_nir
    return trp, rrp, trn, rrn
end

"Per-band result of the roof + screen stack, relative to the same roof with the screens retracted."
@inline function _band(tr, rr, tes, res, tbl, rbl, ues, ubl, bl_on_top)
    es = layer_eff(tes, res, ues); bl = layer_eff(tbl, rbl, ubl)
    up, lo = bl_on_top ? (bl, es) : (es, bl)
    T, R, Ar, A1, A2 = stack3(tr, rr, up[1], up[2], lo[1], lo[2])
    Aes, Abl = bl_on_top ? (A2, A1) : (A1, A2)
    return (S = T / tr, T = T, dR = R - rr, dAr = Ar - (1.0 - tr - rr), Aes = Aes, Abl = Abl)
end

"""
    screen_factors(Ues, Ubl, p) -> NamedTuple

Optics per unit incident global radiation (independent of I2). Ues = energy-screen closure, Ubl = blackout-screen
closure, both fractions 0-1 (the logged EnScr, BlackScr are percentages: divide by 100 before calling).
"""
@inline function screen_factors(Ues, Ubl, p)
    ues = clamp(Ues, 0.0, 1.0); ubl = clamp(Ubl, 0.0, 1.0)
    trp, rrp, trn, rrn = roof_props(p)
    top = p.scr_order > 0.5                         # 1 = blackout screen above the energy screen
    P = _band(trp, rrp, p.scr_es_tau_par, p.scr_es_rho_par, p.scr_bl_tau_par, p.scr_bl_rho_par, ues, ubl, top)
    N = _band(trn, rrn, p.scr_es_tau_nir, p.scr_es_rho_nir, p.scr_bl_tau_nir, p.scr_bl_rho_nir, ues, ubl, top)
    return (P = P, N = N)
end

"""
    screen_solar(I2, Ues, Ubl, p) -> NamedTuple (W m-2 floor)

Spar, Snir: factors on every existing below-screen solar term of that band (1 when retracted).
Qscr: solar absorbed by the screens; Qair = scr_abs_air * Qscr (to the crop air);
Qcov = (1 - scr_abs_air) * Qscr + extra roof absorption (to the cover/structure node); Qrefl: extra reflection out.
"""
@inline function screen_solar(I2, Ues, Ubl, p)
    f = screen_factors(Ues, Ubl, p)
    Ip = p.eta2 * I2; In = p.eta3 * I2
    Qscr = Ip * (f.P.Aes + f.P.Abl) + In * (f.N.Aes + f.N.Abl)
    Qrf  = Ip * f.P.dAr + In * f.N.dAr
    Qair = p.scr_abs_air * Qscr
    return (Spar = f.P.S, Snir = f.N.S, Qscr = Qscr, Qair = Qair, Qcov = Qscr - Qair + Qrf, Qroof = Qrf,
            Qrefl = Ip * f.P.dR + In * f.N.dR,
            Qes = Ip * f.P.Aes + In * f.N.Aes, Qbl = Ip * f.P.Abl + In * f.N.Abl)
end

# =============================================================================
# Note 103 (flag solar_cons, default off): energy-conserving solar column, roof -> screens -> canopy -> floor
# =============================================================================
# Per band (PAR, solar NIR) five elements, all repeated reflections between them (adding method), unit flux from the
# sky:  roof (tau_r, rho_r) | upper screen | lower screen | canopy | opaque floor (reflectance rho_f).
# Canopy = one homogeneous layer over the whole floor (LAI is per m2 floor; no row/aisle split), as Vanthoor (2011)
# eqs 8.31-8.32:  t = exp(-k LAI),  r = rho_inf (1 - t),  a = 1 - t - r  (same from above and below).
#   PAR: k = beta1 (0.7), rho_inf = rho1 (0.07), floor rho2 (0.65)   -- the V2.2 r2/r3 first-pass terms, made exact
#   NIR: k = beta3 (0.27), rho_inf = rho_can_nir (0.35, Vanthoor's canopy NIR REFLECTION), floor 1 - alpha7 (0.5)
# Every absorbed or reflected watt is booked once: canopy, floor (+ mass share), roof (structure share eta1 -> air,
# glass -> cover node), screens (scr_abs_air -> air, rest -> cover node), and reflection out to the sky.
#   canopy + floor + mass + air + cover + out = (eta2 + eta3) I2 = I2.
# Upward light from canopy and floor is traced back through the screens and roof (it is not assumed to escape).
# Lamps (below the screens) and thermal LWIR are outside this calculation.

"""
    column4(t1, r1, t2, r2, t3, r3, t4, r4, rf) -> NamedTuple

Unit flux from above on four symmetric layers over an opaque floor (reflectance rf), repeated reflections included.
Rout + A1 + A2 + A3 + A4 + Af = 1. Dcan = downward flux incident on the top of layer 4 (the canopy).
"""
@inline function column4(t1, r1, t2, r2, t3, r3, t4, r4, rf)
    R4 = r4 + t4 * t4 * rf / (1.0 - r4 * rf)          # reflectance of canopy + floor, seen from above
    R3 = r3 + t3 * t3 * R4 / (1.0 - r3 * R4)
    R2 = r2 + t2 * t2 * R3 / (1.0 - r2 * R3)
    R1 = r1 + t1 * t1 * R2 / (1.0 - r1 * R2)          # whole greenhouse, seen from the sky
    D2 = t1 / (1.0 - r1 * R2)                         # downward fluxes between the layers (D1 = 1)
    D3 = t2 * D2 / (1.0 - r2 * R3)
    D4 = t3 * D3 / (1.0 - r3 * R4)
    D5 = t4 * D4 / (1.0 - r4 * rf)
    return (Rout = R1,
            A1 = (1.0 - t1 - r1) * (1.0 + R2 * D2),   # each layer absorbs from both faces
            A2 = (1.0 - t2 - r2) * (D2 + R3 * D3),
            A3 = (1.0 - t3 - r3) * (D3 + R4 * D4),
            A4 = (1.0 - t4 - r4) * (D4 + rf * D5),
            Af = (1.0 - rf) * D5, Dcan = D4)
end

"Canopy as a homogeneous layer (Vanthoor 2011, eqs 8.31-8.32): (t, r) for extinction k and high-LAI reflection rho_inf."
@inline canopy_layer(k, rho_inf, LAI) = (t = exp(-k * max(LAI, 0.0)); (t, rho_inf * (1.0 - t)))

"""
    solar_column(I2, Ues, Ubl, LAI, p) -> NamedTuple (W m-2 floor)

Ues, Ubl: energy- and blackout-screen closure, fractions 0-1. LAI per m2 floor (the crop model's output).
Qcan  canopy (both bands) | Qflr floor level (split floor / internal mass where the flags are known)
Qair  crop air: roof structure share + scr_abs_air of the screens | Qcov cover/structure node: roof glass + rest of
the screens | Qout reflected to the sky | ParTop downward sun PAR on the canopy top (W m-2 PAR; plant light).
"""
@inline function solar_column(I2, Ues, Ubl, LAI, p)
    ues = clamp(Ues, 0.0, 1.0); ubl = clamp(Ubl, 0.0, 1.0)
    trp, rrp, trn, rrn = roof_props(p)
    top = p.scr_order > 0.5
    esP = layer_eff(p.scr_es_tau_par, p.scr_es_rho_par, ues); blP = layer_eff(p.scr_bl_tau_par, p.scr_bl_rho_par, ubl)
    esN = layer_eff(p.scr_es_tau_nir, p.scr_es_rho_nir, ues); blN = layer_eff(p.scr_bl_tau_nir, p.scr_bl_rho_nir, ubl)
    uP, lP = top ? (blP, esP) : (esP, blP)
    uN, lN = top ? (blN, esN) : (esN, blN)
    cP = canopy_layer(p.beta1, p.rho1, LAI)
    cN = canopy_layer(p.beta3, p.rho_can_nir, LAI)
    P = column4(trp, rrp, uP[1], uP[2], lP[1], lP[2], cP[1], cP[2], p.rho2)
    N = column4(trn, rrn, uN[1], uN[2], lN[1], lN[2], cN[1], cN[2], 1.0 - p.alpha7)
    Ip = p.eta2 * I2; In = p.eta3 * I2
    arP = 1.0 - trp - rrp; arN = 1.0 - trn - rrn
    fsP = arP > 0.0 ? min(p.eta1 / arP, 1.0) : 0.0      # structure share of the roof absorption (-> air)
    fsN = arN > 0.0 ? min(p.eta1 / arN, 1.0) : 0.0
    Qroof = Ip * P.A1 + In * N.A1
    Qst   = Ip * P.A1 * fsP + In * N.A1 * fsN
    Qes   = top ? Ip * P.A3 + In * N.A3 : Ip * P.A2 + In * N.A2
    Qbl   = top ? Ip * P.A2 + In * N.A2 : Ip * P.A3 + In * N.A3
    Qscr  = Qes + Qbl
    Qsa   = p.scr_abs_air * Qscr
    return (Qcan = Ip * P.A4 + In * N.A4, QcanP = Ip * P.A4, QcanN = In * N.A4,
            Qflr = Ip * P.Af + In * N.Af, QflrP = Ip * P.Af, QflrN = In * N.Af,
            Qair = Qst + Qsa, Qst = Qst, Qsa = Qsa, Qcov = (Qroof - Qst) + (Qscr - Qsa), Qroof = Qroof,
            Qscr = Qscr, Qes = Qes, Qbl = Qbl, Qout = Ip * P.Rout + In * N.Rout,
            ParTop = Ip * P.Dcan, NirTop = In * N.Dcan)
end
