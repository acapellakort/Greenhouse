# =============================================================================
# solar.jl -- energy-conserving solar column (V2.2 Notes 101-105, flags solar_cons + par_below)
# =============================================================================
# Per band (PAR = eta2 Iglob, solar NIR = eta3 Iglob), five elements from the sky down, all repeated reflections
# (adding method, Vanthoor 2011 eqs 8.14-8.17):
#     roof (tau_r, rho_r) | energy screen | blackout screen | canopy | opaque floor (reflectance rho_f)
# Screen deployed by the fraction u:      tau = 1 - u (1 - tau_s),  rho = u rho_s          (Vanthoor 8.14-8.17)
# Canopy (one homogeneous layer):         t = exp(-k LAI),  r = rho_inf (1 - t),  a = 1 - t - r   (Vanthoor 8.31-8.32)
#     PAR: k = beta1, rho_inf = rho1, floor rho2;   NIR: k = beta3, rho_inf = rho_can_nir, floor 1 - alpha7
# Roof: tau_r = (1 - eta1) tau1 for both bands; rho_r = 1 - tau_r - eta1 - cov_abs, so that the roof absorbs eta1
# (structure, booked to the air) + cov_abs (glass, booked to the cover node).
# Booking (W m-2 floor):  canopy Qcan | floor level Qflr (floor (1 - f_mass), internal mass f_mass) | air Qair =
# roof structure + scr_abs_air of the screens | cover Qcov = roof glass + rest of the screens | sky Qout.
#     Qcan + Qflr + Qair + Qcov + Qout = Iglob   (eta2 + eta3 = 1)
# Lamps hang below the screens and are not part of this calculation; long-wave radiation is not handled here.

"Solar fluxes of one time instant [W m-2 floor] (all fields Float64; isbits, so rhs! does not allocate)."
struct SolarFluxes
    Qcan::Float64     # canopy, both bands
    QcanP::Float64    # canopy, PAR
    QcanN::Float64    # canopy, NIR
    Qflr::Float64     # floor level (floor + internal mass), both bands
    QflrP::Float64
    QflrN::Float64
    Qair::Float64     # crop air: roof structure + screen share
    Qst::Float64      # roof structure part of Qair
    Qsa::Float64      # screen part of Qair
    Qcov::Float64     # cover node: roof glass + rest of the screen absorption
    Qroof::Float64    # total roof absorption
    Qscr::Float64     # total screen absorption
    Qes::Float64      # energy screen absorption
    Qbl::Float64      # blackout screen absorption
    Qout::Float64     # reflected back to the sky
    ParTop::Float64   # downward PAR on the canopy top [W m-2 PAR]
    NirTop::Float64   # downward NIR on the canopy top
    ParBot::Float64   # PAR reaching the canopy underside (floor-reflected, after all returns)
end

"Screen deployed by the fraction u: effective (tau, rho)."
@inline layer_eff(tau, rho, u) = (1.0 - u * (1.0 - tau), u * rho)

"Roof optics (tau_par, rho_par, tau_nir, rho_nir); NIR equal to PAR."
@inline function roof_props(p)
    trp = (1.0 - p.eta1) * p.tau1
    rrp = max(1.0 - trp - p.eta1 - p.cov_abs, 0.0)
    return trp, rrp, trp, rrp
end

"""
    column4(t1, r1, t2, r2, t3, r3, t4, r4, rf) -> NamedTuple

Unit flux from above on four symmetric layers over an opaque floor (reflectance rf), all repeated reflections.
Rout + A1 + A2 + A3 + A4 + Af = 1. Dcan = downward flux on the top of layer 4, Dflr = downward flux on the floor.
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
            Af = (1.0 - rf) * D5, Dcan = D4, Dflr = D5)
end

"Canopy as a homogeneous layer: (t, r) for extinction k and high-LAI reflection rho_inf."
@inline canopy_layer(k, rho_inf, LAI) = (t = exp(-k * max(LAI, 0.0)); (t, rho_inf * (1.0 - t)))

"""
    solar_column(Iglob, Ues, Ubl, LAI, p) -> SolarFluxes

Ues, Ubl: energy- and blackout-screen closure (fractions, clamped to 0-1); the energy screen is the upper screen.
"""
@inline function solar_column(Iglob, Ues, Ubl, LAI, p)
    ues = clamp(Ues, 0.0, 1.0); ubl = clamp(Ubl, 0.0, 1.0)
    trp, rrp, trn, rrn = roof_props(p)
    esP = layer_eff(p.scr_es_tau_par, p.scr_es_rho_par, ues); blP = layer_eff(p.scr_bl_tau_par, p.scr_bl_rho_par, ubl)
    esN = layer_eff(p.scr_es_tau_nir, p.scr_es_rho_nir, ues); blN = layer_eff(p.scr_bl_tau_nir, p.scr_bl_rho_nir, ubl)
    cP = canopy_layer(p.beta1, p.rho1, LAI)
    cN = canopy_layer(p.beta3, p.rho_can_nir, LAI)
    P = column4(trp, rrp, esP[1], esP[2], blP[1], blP[2], cP[1], cP[2], p.rho2)
    N = column4(trn, rrn, esN[1], esN[2], blN[1], blN[2], cN[1], cN[2], 1.0 - p.alpha7)
    Ip = p.eta2 * Iglob; In = p.eta3 * Iglob
    arP = 1.0 - trp - rrp; arN = 1.0 - trn - rrn
    fsP = arP > 0.0 ? min(p.eta1 / arP, 1.0) : 0.0      # structure share of the roof absorption (-> air)
    fsN = arN > 0.0 ? min(p.eta1 / arN, 1.0) : 0.0
    Qroof = Ip * P.A1 + In * N.A1
    Qst   = Ip * P.A1 * fsP + In * N.A1 * fsN
    Qes   = Ip * P.A2 + In * N.A2
    Qbl   = Ip * P.A3 + In * N.A3
    Qscr  = Qes + Qbl
    Qsa   = p.scr_abs_air * Qscr
    return SolarFluxes(Ip * P.A4 + In * N.A4, Ip * P.A4, In * N.A4, Ip * P.Af + In * N.Af, Ip * P.Af, In * N.Af,
                       Qst + Qsa, Qst, Qsa, (Qroof - Qst) + (Qscr - Qsa), Qroof, Qscr, Qes, Qbl,
                       Ip * P.Rout + In * N.Rout, Ip * P.Dcan, In * N.Dcan, Ip * p.rho2 * P.Dflr)
end

"""
    par_both_faces_factor(LAI, p) -> g >= 1

Crop path (measured canopy-top PAR): PAR on both canopy faces per unit PAR transmitted by the roof, from the PAR
column with the screens retracted, g = (D_top + rho2 D_floor) / tau_roof.
"""
@inline function par_both_faces_factor(LAI, p)
    trp, rrp, _, _ = roof_props(p)
    cP = canopy_layer(p.beta1, p.rho1, LAI)
    c  = column4(trp, rrp, 1.0, 0.0, 1.0, 0.0, cP[1], cP[2], p.rho2)
    return (c.Dcan + p.rho2 * c.Dflr) / trp
end
