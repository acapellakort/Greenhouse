# =============================================================================
# photosynthesis.jl  --  FvCB (Farquhar-von Caemmerer-Berry) C3 model
# =============================================================================
#
# Faithful port of 2.0's `AcropFast`. Changes are mechanical only:
#   * crop constants are destructured from the concrete NamedTuple `p`;
#   * `minimum([Ar_j, Ar_c, Ar_p])` -> `min(Ar_j, Ar_c, Ar_p)` to avoid a heap
#     allocation on every call (numerically identical);
#   * stomatal conductance kept as the constant effective conductance `g_eff`
#     (Bush 2023) exactly as in 2.0.
# Result is bit-for-bit equal to 2.0 for the same inputs.

"Effective (stomatal+mesophyll) conductance [mol m^-2 s^-1] -- Bush (2023)."
@inline g_eff() = (0.3 * 0.25) / (0.3 + 0.25)

"""
    jarvis_gs(Iw, C, VPD, T, p) -> gs   [mol m^-2 s^-1]

Jarvis-type multiplicative STOMATAL conductance: gs_max scaled by independent
light, VPD and CO2 stress factors, floored at gs_min (cuticular). Combined in
series with a fixed mesophyll conductance by `stomatal_conductance`.
`Iw` uses the same light units as `assimilation` (incident W m^-2); `C` is CO2
[mg m^-3]; `VPD` [kPa]. Temperature response is left to the FvCB biochemistry
(Vcmax/Jmax) and not double-counted here.
"""
@inline function jarvis_gs(Iw, C, VPD, T, p)
    (; gs_max, gs_I_half, gs_D0, gs_Ca_ref, gs_Ca_scale, gs_min) = p
    ppm = C / 1.83                                   # mg m^-3 -> ppm (~20C)
    fI  = Iw / (Iw + gs_I_half)                      # light saturation
    fD  = 1.0 / (1.0 + max(VPD, 0.0) / gs_D0)        # closes under high VPD
    fC  = 1.0 / (1.0 + max(ppm - gs_Ca_ref, 0.0) / gs_Ca_scale)  # closes under CO2 enrichment
    return max(gs_min, gs_max * fI * fD * fC)
end

"Total conductance: constant g_eff (stomata_model=0), or Jarvis stomatal in series with fixed mesophyll (=1)."
@inline function stomatal_conductance(Iw, C, VPD, T, p)
    p.stomata_model == 0.0 && return g_eff()
    gs = jarvis_gs(Iw, C, VPD, T, p)          # stomatal component (Jarvis-modulated)
    gm = p.gs_mesophyll                        # fixed mesophyll conductance
    return (gs * gm) / (gs + gm)               # total = stomatal in series with mesophyll
end

"""
    assimilation(T, C, Iw, LAI, p) -> A   [mg CO2 m^-2 s^-1]

Canopy assimilation from the FvCB model: the minimum of the RuBisCO-limited,
electron-transport(light)-limited and product-limited rates, scaled by LAI.
`T` [K], `C` CO2 [mg m^-3], `Iw` incident radiation [W m^-2], `p` the crop
parameter NamedTuple.

NOTE (fidelity): 2.0 calls this with `LAI = I9` (radiation) rather than the
true leaf-area index -- a known quirk kept here so V2.2 reproduces 2.0 exactly.
The forthcoming growth model should pass the real LAI instead; do that by
changing the `assimilation(...)` call in forward_map.jl / rhs!.
"""
@inline function assimilation(T, C, Iw, LAI, p; VPD = 0.0)
    (; O_a, Sco25, E_Soc, Rgas, alpha, J_max, theta,
       V_cmax25, Q10_Vcmax, Rd_day, K_C25, Q10_KC, K_O25, Q10_KO) = p

    I   = 4.6e-6 * Iw          # W m^-2 -> mol photons m^-2 s^-1
    C_i = 0.509e-6 * C         # mg m^-3 -> mol m^-3

    # CO2 compensation point [mol m^-3]
    Gamma_st = (0.5 * O_a) / (Sco25 * exp(((T - 298.15) / T) * (E_Soc) / (298.15 * Rgas)))

    # electron transport rate [mol m^-2 s^-1]
    J_t = ((alpha * I + J_max) - sqrt((alpha * I + J_max)^2 - 4 * theta * alpha * I * J_max)) / (2 * theta)

    g = stomatal_conductance(Iw, C, VPD, T, p)

    # temperature-adjusted Vcmax and Michaelis-Menten constants
    V_cmax = V_cmax25 * Q10_Vcmax^(0.1 * (T - 298.15)) / (1 + exp(0.128 * (T - 315.15)))
    K_C = K_C25 * Q10_KC^(0.1 * (T - 298.15))
    K_O = K_O25 * Q10_KO^(0.1 * (T - 298.15))

    # RuBisCO-limited
    pc = -(V_cmax + g * (C_i + K_C * (1 + O_a / K_O)) - Rd_day)
    qc = g * (V_cmax * max(C_i - Gamma_st, 0) - (C_i + K_C * (1 + O_a / K_O)) * Rd_day)
    disc_c = pc^2 - 4 * qc
    Ar_c = disc_c >= 0 ? 0.5 * (-pc - sqrt(disc_c)) : -Rd_day

    # light-limited
    pj = -0.25 * J_t + Rd_day - g * (C_i + 2.0 * Gamma_st)
    qj = 0.25 * g * max(C_i - Gamma_st, 0) * J_t - g * (C_i + 2 * Gamma_st) * Rd_day
    disc_j = pj^2 - 4 * qj
    Ar_j = disc_j >= 0 ? 0.5 * (-pj - sqrt(disc_j)) : -Rd_day

    # product-limited
    Ar_p = 1.5 * V_cmax - Rd_day

    # mol m^-2 s^-1 -> mg m^-2 s^-1  (molar mass CO2 = 44.010 g/mol)
    return 44010 * LAI * min(Ar_j, Ar_c, Ar_p)
end
