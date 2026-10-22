# =============================================================================
# photosynthesis.jl -- FvCB (Farquhar-von Caemmerer-Berry) canopy assimilation
# =============================================================================
# Port of V2.2 `assimilation` with the constant effective CO2 conductance (stomata_model = 0) and the leaf dark
# respiration scaled by Rd_scale (Q10 = 2 on the canopy temperature in the climate model).

"Effective light-intercepting leaf area (Beer's law): (1 - exp(-k LAI)) / k."
@inline lai_eff(lai, k) = (1.0 - exp(-k * lai)) / k

"""
    assimilation(T, C, Iw, Leff, p; Rd_scale = 1.0) -> A   [mg CO2 m-2 s-1]

T [K] leaf (air) temperature, C [mg m-3] CO2, Iw [W m-2] PAR, Leff effective leaf area. Minimum of the
RuBisCO-, light- and product-limited rates times Leff.
"""
@inline function assimilation(T, C, Iw, Leff, p; Rd_scale = 1.0)
    (; O_a, Sco25, E_Soc, Rgas, alpha, J_max, theta,
       V_cmax25, Q10_Vcmax, K_C25, Q10_KC, K_O25, Q10_KO) = p
    Rd_day = p.Rd_day * Rd_scale

    I   = 4.6e-6 * Iw          # W m-2 -> mol photons m-2 s-1
    C_i = 0.509e-6 * C         # mg m-3 -> mol m-3

    Gamma_st = (0.5 * O_a) / (Sco25 * exp(((T - 298.15) / T) * (E_Soc) / (298.15 * Rgas)))   # CO2 compensation point
    J_t = ((alpha * I + J_max) - sqrt((alpha * I + J_max)^2 - 4 * theta * alpha * I * J_max)) / (2 * theta)
    g = (p.g_s * p.g_m) / (p.g_s + p.g_m)              # constant effective conductance (Bush 2023)

    V_cmax = V_cmax25 * Q10_Vcmax^(0.1 * (T - 298.15)) / (1 + exp(0.128 * (T - 315.15)))
    K_C = K_C25 * Q10_KC^(0.1 * (T - 298.15))
    K_O = K_O25 * Q10_KO^(0.1 * (T - 298.15))

    pc = -(V_cmax + g * (C_i + K_C * (1 + O_a / K_O)) - Rd_day)                         # RuBisCO-limited
    qc = g * (V_cmax * max(C_i - Gamma_st, 0) - (C_i + K_C * (1 + O_a / K_O)) * Rd_day)
    disc_c = pc^2 - 4 * qc
    Ar_c = disc_c >= 0 ? 0.5 * (-pc - sqrt(disc_c)) : -Rd_day

    pj = -0.25 * J_t + Rd_day - g * (C_i + 2.0 * Gamma_st)                              # light-limited
    qj = 0.25 * g * max(C_i - Gamma_st, 0) * J_t - g * (C_i + 2 * Gamma_st) * Rd_day
    disc_j = pj^2 - 4 * qj
    Ar_j = disc_j >= 0 ? 0.5 * (-pj - sqrt(disc_j)) : -Rd_day

    Ar_p = 1.5 * V_cmax - Rd_day                                                        # product-limited

    return 44010 * Leff * min(Ar_j, Ar_c, Ar_p)       # mol m-2 s-1 -> mg CO2 m-2 s-1
end
