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
    trp = (1.0 - p.eta1) * p.tau1
    trn = p.roof_tau_nir < 0.0 ? trp : p.roof_tau_nir
    top = p.scr_order > 0.5                         # 1 = blackout screen above the energy screen
    P = _band(trp, p.roof_rho_par, p.scr_es_tau_par, p.scr_es_rho_par, p.scr_bl_tau_par, p.scr_bl_rho_par, ues, ubl, top)
    N = _band(trn, p.roof_rho_nir, p.scr_es_tau_nir, p.scr_es_rho_nir, p.scr_bl_tau_nir, p.scr_bl_rho_nir, ues, ubl, top)
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
