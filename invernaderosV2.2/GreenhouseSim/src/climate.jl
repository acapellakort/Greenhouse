# =============================================================================
# climate.jl  --  greenhouse climate ODE right-hand side (type-stable)
# =============================================================================
#
# Faithful port of 2.0's `rhs_fast` (in climate_functionsV3_NoUnits.jl). The
# ONLY change is mechanical: every physical constant that was previously read
# from a non-`const` global is now destructured from the concrete parameter
# NamedTuple `p` at the top of the function via `(; name, ...) = p`. The math
# is byte-for-byte identical, so results match 2.0 (see test/benchmark.jl for
# the regression check) -- but the function is now fully inferable and
# allocation-free.

# --- psychrometric helpers (unchanged from 2.0) ------------------------------
function Pws(T1)
    return ifelse(T1 > 273.15,
        611.21 * exp((18.678 - ((T1 - 273.15) / 234.5)) * ((T1 - 273.15) / (257.14 + T1 - 273.15))),
        611.15 * exp((23.036 - ((T1 - 273.15) / 333.7)) * ((T1 - 273.15) / (279.82 + T1 - 273.15))))
end

"Relative humidity [%] from air temperature [K] and vapour pressure [Pa]."
rhf(T1, V1) = 100 * V1 / Pws(T1)

"Vapour-pressure deficit [Pa]."
VPDf(VP_sat, rh) = VP_sat * (1 - rh / 100)

"""
    cover_T(T2, I4, I5, I2, p) -> Tc

Algebraic cover temperature (Note 29 cover node). With solar_split = 1 (Note 46) the cover also absorbs
cov_abs * I2 of the global radiation:  Tc = (kR(T2 + I4) + h_cover_out*I5 + R_cov) / (2kR + h_cover_out).
"""
@inline function cover_T(T2, I4, I5, I2, p)
    kR  = 4.0 * p.epsil3 * p.sigma * (0.5 * (T2 + I4))^3
    Rc  = (hasproperty(p, :solar_split) && p.solar_split > 0.5) ? p.cov_abs * I2 : 0.0
    return (kR * (T2 + I4) + p.h_cover_out * I5 + Rc) / (2.0 * kR + p.h_cover_out)
end

"""
    cover_latent(Tbase, T2, V1, U1, kR, p) -> Tc

Note 81 (only with cover_cond = 1 and the cover node on): the latent heat of the vapour condensing on the cover
(Vanthoor MV_AirCov, per m2 cover: 6.4e-9 * 1.7|T2 - Tc|^0.33 * (VP - Psat(Tc)) * share of crop air in contact)
is released IN the cover, which warms it and limits condensation. The cover node is linear in Tc apart from this
term, so Tc = Tbase + gamma2 * MV(Tc) / (2kR + h_cover_out), solved by fixed point (contraction ~0.1, 4 passes).
The vapour removed from the air is the same MV (climate_rhs, p10), so water and energy stay consistent.
"""
@inline function cover_latent(Tbase, T2, V1, U1, kR, p)
    D    = 2.0 * kR + p.h_cover_out
    fcov = 1.0 - U1 * (1.0 - p.vent_scr_perm)
    Tc   = Tbase
    for _ in 1:4
        dv = V1 - Pws(Tc)
        Tc = Tbase + (dv > 0.0 ? p.gamma2 * fcov * 6.4e-9 * 1.7 * abs(T2 - Tc)^0.33 * dv / D : 0.0)
    end
    return Tc
end
cond_on(p) = hasproperty(p, :cover_cond) && p.cover_cond > 0.5

"""
    floor_full_fluxes(T1, T2, T5, I1, I2, I4, I5, U1, p) -> (dT1_corr, dT2_corr, Q_floor)

Note 46. Re-books the floor terms of climate_rhs so the floor energy balance closes (floor_full = 1):
  r11 (floor FIR to sky)  was taken from the AIR       -> given back to the air, taken from the floor
  r14 (canopy FIR to floor) used T2 and was taken from the canopy only -> recomputed with T1, given to the floor
  h11 (air-floor conduction) was taken from the AIR only -> given to the floor
With solar_split = 1 the floor also absorbs the transmitted sun not intercepted by the canopy.
Returns corrections to dT1, dT2 [K/s] and the extra net flux into the floor [W m^-2].
"""
@inline function floor_full_fluxes(T1, T2, T5, I1, I2, I4, I5, U1, p; V1 = NaN)
    (; epsil2, epsil3, epsil4, sigma, gamma1, phi1, tau1, tau3, beta1, beta3, eta1, eta2, eta3,
       rho2, alpha1, alpha5, alpha7, rho3, phi2, nu7, nu8, use_cover_node) = p
    Qf = 0.0; c1 = 0.0; c2 = 0.0
    if hasproperty(p, :floor_full) && p.floor_full > 0.5
        Tc  = use_cover_node > 0.5 ? cover_T(T2, I4, I5, I2, p) : I4
        if use_cover_node > 0.5 && cond_on(p) && isfinite(V1)      # Note 81: same cover T as climate_rhs
            Tc = cover_latent(Tc, T2, V1, U1, 4.0 * epsil3 * sigma * (0.5 * (T2 + I4))^3, p)
        end
        g2  = 1 - U1 * (1 - tau3)
        g3  = tau1 * g2 * (1 - 0.49 * pi * gamma1 * phi1) * exp(-beta3 * I1)
        fcf = epsil4 * epsil2 * (1 - 0.49 * pi * gamma1 * phi1) * sigma
        r11  = epsil4 * g3 * epsil3 * sigma * (T5^4 - Tc^4)
        r14o = fcf * (T2^4 - T5^4)
        r14n = fcf * (T1^4 - T5^4)
        h11  = 2 * nu7 * (T2 - T5) / (phi2 + nu8)
        c2 = r11 / (phi2 * rho3 * alpha5)
        c1 = (r14o - r14n) / (alpha1 * I1)
        Qf = h11 + r14n - r11
    end
    if hasproperty(p, :solar_split) && p.solar_split > 0.5
        Qf += (1 - eta1) * tau1 * I2 * (eta2 * exp(-beta1 * I1) * (1 - rho2) + eta3 * exp(-beta3 * I1) * alpha7)
    end
    return c1, c2, Qf
end

"""
    climate_rhs(T1,T2,C1,V1, I1..I11, A, U1..U10,U12, p) -> (dT1, dT2, dC1, dV1)

Right-hand sides of the four greenhouse-climate ODEs (canopy temperature,
air temperature, CO2, vapour pressure). `p` is the concrete parameter
NamedTuple from `load_params`. `A` is the canopy assimilation rate.
"""
@inline function climate_rhs(T1, T2, C1, V1,
        I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A,
        U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12, p)

    # --- all physical constants, read once from the concrete NamedTuple ------
    (; alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9, alpha12,
       beta1, beta2, beta3,
       delta1, delta2, delta3, delta4, delta5, delta6, delta7,
       epsil1, epsil2, epsil3, epsil4, epsil5, epsil6, eps_screen,
       eta1, eta2, eta3, eta5, eta6, eta7, eta8, eta10, eta11, eta12, eta13, eta14, eta15, eta16,
       gamma, gamma1, gamma2, gamma3, gamma4, gamma5,
       lamb1, lamb2, lamb4,
       n_pipes, nu1, nu2, nu3, nu4, nu4CO2, nu5, nu6, nu7, nu8,
       omega1, omega2,
       phi1, phi2, phi5, phi6, phi7, phi8, phi9,
       psi1, psi2,
       rho1, rho2, rho3,
       sigma,
       tau1, tau2, tau3,
       use_cover_node, h_cover_out,
       nu4_min,
       U8_min_leak) = p

    # --- global quantities ---------------------------------------------------
    T_ThmScr = 0.5 * (I5 + T2)
    VPsat_ThmScr = Pws(T_ThmScr)
    VPsat_Air = Pws(T2)
    q2 = Pws(T1)                        # VP of canopy

    # stomatal-resistance auxiliary functions
    q7 = (1 + exp(1 / gamma5 * (I9 - delta1)))^-1
    q8 = delta4 + (delta5 - delta4) * q7
    q9 = delta6 + (delta7 - delta6) * q7
    q10 = (I9 + delta2) / (I9 + delta3)

    q6 = Pws(I6)                        # VP of mechanical cooling

    q5 = 1 + q9 * ((q2 - V1)^2)         # levels large VP deviations
    q4 = 1 + q8 * ((C1 - 200)^2)        # levels large CO2 values
    # --- stomata_fix (default-off: identical to V2.2) --------------------------
    # q7 is ~1 by DAY (gamma5 < 0), but q8/q9 above give weight q7 to the NIGHT
    # coefficients (delta5, delta7) -> day/night swapped w.r.t. Vanthoor (2011).
    # Also delta4/delta5 are per ppm^2 while C1 is mg/m3 (~1.8 x ppm).
    # Fix: day coefficients weighted by q7, and C1 converted to ppm (mg/m3 * T / 536.4).
    if hasproperty(p, :stomata_fix) && p.stomata_fix > 0.5
        q8f = delta5 + (delta4 - delta5) * q7
        q9f = delta7 + (delta6 - delta7) * q7
        q5  = 1 + q9f * ((q2 - V1)^2)
        q4  = 1 + q8f * ((C1 * T2 / 536.4 - 200)^2)
    end
    q3 = gamma4 * q10 * q4 * q5         # canopy stomatal resistance
    q1 = (2 * rho3 * alpha5 * I1) / (gamma * gamma2 * (gamma3 + q3))  # VEC Can-Air

    # ventilation rates
    n1 = nu1 * (1 - eta10 * U5)
    n2 = nu3 * U6
    n3 = nu2 * (1 - eta11 * U5)
    f5 = n1 * n2 * I8 * sqrt(n3) / (2 * alpha6)

    f6 = max(I8 < 0.25 ? 0.25 * nu4 : nu4 * I8, nu4_min)
    f6CO2 = max(I8 < 0.25 ? 0.25 * nu4CO2 : nu4CO2 * I8, nu4_min)

    f7t = max(omega1 * nu6 * (T2 - I5) / (T2 + I5) + n3 * I8^2, 0.0)
    # Ventilation RATE is a non-negative air-exchange rate (buoyancy magnitude);
    # the DIRECTION of heat/CO2/vapour transport is already carried by (T2-I5),
    # (C1-I10), (V1/T2-I11/I5) in the balances. A previous `* sign(T2-I5)` factor
    # made the rate negative when inside was colder than outside, turning the CO2
    # and vapour ventilation terms into a positive feedback (blow-up) once the
    # leakage is set to a realistic value. Removed. Identical whenever T2 >= I5
    # (all venting in the summer calibration), so calibrated results are unchanged.
    # --- roof-vent effectiveness flags (default-off: identical to V2.2) --------
    # vent_exp:    effective opening = U8^vent_exp (real vent flow vs opening % is concave;
    #              vent_exp = 1 -> linear, skipped by an explicit branch so it is exact).
    # vent_screen: a closed energy screen (U1) throttles roof ventilation of the CROP air,
    #              factor (1 - U1*(1 - vent_scr_perm)): the vents then mostly ventilate the
    #              space above the screen. One-node proxy for Vanthoor's top compartment.
    U8e = max(U8, U8_min_leak)
    if hasproperty(p, :vent_exp) && p.vent_exp != 1.0
        U8e = U8e^p.vent_exp
    end
    if hasproperty(p, :vent_screen) && p.vent_screen > 0.5
        U8e = U8e * (1.0 - U1 * (1.0 - p.vent_scr_perm))
    end
    f7 = (U8e * nu5 * n1) / (2.0 * alpha6) * sqrt(f7t)

    if eta7 >= eta8
        f2 = eta6 * f5 + 0 * 0.5 * f6
    else
        f2 = eta6 * (U1 * f5) + 0 * 0.5 * f6
    end
    f3 = U7 * phi8 / alpha6
    f4 = eta6 * f7 + 0 * 0.5 * f6

    # PAR / NIR / FIR radiation
    g1 = 0.49 * (1 - exp(-beta3 * I1))
    r4 = (1 - eta1) * tau1 * eta2 * I2
    r2 = r4 * (1 - rho1) * (1 - exp(-beta1 * I1))
    r3 = r4 * exp(-beta1 * I1) * rho2 * (1 - rho1) * (1 - exp(-beta2 * I1))

    g2 = 1 - U1 * (1 - tau3)
    g3 = tau1 * g2 * (1 - 0.49 * pi * gamma1 * phi1) * exp(-beta3 * I1)

    # --- algebraic cover node (default-off: use_cover_node=0 → I4_eff = I4) ----
    # Quasi-static cover energy balance (linearised around T_ref = (T2+I4)/2):
    #   kR*(T2-Tc) = kR*(Tc-I4) + h_cover_out*(Tc-I5)
    #   → Tc = (kR*(T2+I4) + h_cover_out*I5) / (2kR + h_cover_out)
    # When use_cover_node=0: I4_eff = I4 exactly (bit-for-bit with V2.2).
    T_ref_cov = 0.5 * (T2 + I4)
    kR_cov    = 4.0 * epsil3 * sigma * T_ref_cov^3
    I4_eff    = use_cover_node > 0.5 ?
        (kR_cov * (T2 + I4) + h_cover_out * I5) / (2.0 * kR_cov + h_cover_out) : I4
    # Note 46 (default-off): solar_split = 1 -> the global radiation absorbed by the cover (cov_abs * I2) heats
    # the cover node instead of the greenhouse air (see r8 below). Only active with use_cover_node = 1.
    if hasproperty(p, :solar_split) && p.solar_split > 0.5 && use_cover_node > 0.5
        I4_eff = cover_T(T2, I4, I5, I2, p)
    end
    # Note 81: with cover condensation on, its latent heat warms the cover (self-limiting condensation)
    if use_cover_node > 0.5 && cond_on(p)
        I4_eff = cover_latent(I4_eff, T2, V1, U1, kR_cov, p)
    end

    SkyFIRBoltzmann = epsil3 * sigma * (T2^4 - I4_eff^4)
    r11 = epsil4 * g3 * epsil3 * sigma * (I7^4 - I4_eff^4)
    r12 = eps_screen * (U1 * tau2) * SkyFIRBoltzmann   # energy screen (low-e); was epsil5=1 (full-emissivity)
    r13 = epsil6 * SkyFIRBoltzmann
    # Note 79 (default 1 = V2.2): exterior cover area per floor area (AGC compartment: roof + one outside wall ~1.6).
    # The air-cover exchange r13 scales with the cover area; the algebraic cover temperature does not (both of its
    # sides scale alike).
    cov_area = hasproperty(p, :cover_area) ? p.cover_area : 1.0
    r13 = cov_area * r13
    r14 = epsil4 * epsil2 * (1 - 0.49 * pi * gamma1 * phi1) * sigma * (T2^4 - I7^4)
    r15 = epsil2 * epsil5 * (1 - exp(-beta2 * I1)) * g2 * sigma * (T1^4 - I4_eff^4)

    # --- dT1: canopy temperature ---------------------------------------------
    r1 = r2 + r3
    r5 = (1 - eta1) * alpha2 * (eta3 * I2 + eta14 * alpha12 * U12)
    r6 = alpha3 * epsil1 * epsil2 * g1 * sigma * (I3^4 - T1^4)
    h1 = 2 * alpha4 * I1 * (T1 - T2)
    l1 = gamma2 * q1 * (q2 - V1)
    r7 = (g1 / 0.49) * epsil2 * epsil3 * g2 * sigma * (T1^4 - I4_eff^4)

    dT1 = (1.0 / (alpha1 * I1)) * (r1 + r5 + r6 - h1 - l1 - r7 - r14)

    # --- dT2: air temperature ------------------------------------------------
    h2 = (U2 * phi7 / alpha6) * (rho3 * alpha5 * I5 - gamma2 * rho3 * (eta5 * (phi5 - phi6)))
    h3 = (U3 * lamb1 * lamb2 / alpha6) * (I6 - T2) / (T2 - I6 + (6.4e-9 * gamma2 * (V1 - q6)))  # (kept for fidelity; unused below)
    h4 = n_pipes * (1.99 * pi * phi1 * gamma1 * abs((I3 - T2))^0.32) * (I3 - T2)
    # h7: sensible-heat loss by ventilation. Original 2.0 (and the first V2.2 port)
    # omitted f4 (roof/buoyancy vent) here, though the CO2 (o5) and vapour (p5)
    # balances both include it. Consequence: on calm days the wind-driven side vent
    # f2 -> 0 and the greenhouse could not shed heat -> runaway to ~40C. A ventilation
    # air-exchange rate carries sensible heat regardless of what drives it, so the roof
    # vent must appear here too. Added f4. (Changes the venting-day forward map, so
    # climate inference should be re-checked -- unlike the sign fix, this is NOT identical.)
    h7 = rho3 * alpha5 * (f2 + f4 + f3 + 0.5 * f6) * (T2 - I5)
    h11 = 2 * nu7 * (T2 - I7) / (phi2 + nu8)
    h12 = (eta15 + eta16) * alpha12 * U12
    r8 = I2 * (eta1 * (tau1 * eta2 + (alpha2 + alpha7) * eta3) + (alpha8 * eta2 + alpha9 * eta3))
    # Note 46 (default-off): V2.0 adds (alpha8*eta2 + alpha9*eta3)*I2 = 1.0*I2 ("cover absorption", values not
    # given -> 1) to the AIR, so the air receives ~108 % of the outside global radiation on top of the canopy share.
    # In Vanthoor (2011) the cover-absorbed part heats the cover (R_Glob_Cov), not the air. solar_split = 1 keeps
    # only R_Glob_Air here; the cover part goes to the cover node (cover_T) and a floor share to the floor (forward_map).
    if hasproperty(p, :solar_split) && p.solar_split > 0.5
        r8 = I2 * eta1 * (tau1 * eta2 + (alpha2 + alpha7) * eta3)
    end
    r10 = r11 + r12 + r13
    l2 = gamma2 * U9 * phi9 / alpha6

    dT2 = (1 / (phi2 * rho3 * alpha5)) * (h1 + h2 + h4 - h7 - h11 + h12 + r8 - r10 - l2)

    # --- dV1: vapour pressure ------------------------------------------------
    p1 = q1 * (q2 - V1)
    p2 = rho3 * (U2 * phi7 / alpha6) * (eta5 * (phi5 - phi6) + phi6)
    p3 = U9 * phi9 / alpha6
    p4 = eta12 * U4 * lamb4 / alpha6
    p5 = (psi1 / omega2) * (V1 / T2 - I11 / I5) * (f2 + f3 + f4 + f6)
    p6 = (U2 * phi7 / alpha6) * (psi1 / omega2) * (V1 / T2)
    p7 = V1 < q6 ? 0.0 : 6.4e-9 * (V1 - q6)
    p8 = V1 < VPsat_ThmScr ? 0.0 : 6.4e-9 * 1.7 * (abs(T2 - VPsat_ThmScr))^0.33 * (V1 - VPsat_ThmScr)
    p9 = V1 < VPsat_Air ? 0.0 : 6.4e-9 * 1.2 * 1000 * ((abs(T1 - T2))^0.33) * (V1 - VPsat_Air)

    # --- vapour_fix (default-off: identical to V2.2) ---------------------------
    # Vanthoor (2011): MV = 6.4e-9 * HEC * (VP_air - VP_surface), HEC in W m^-2 K^-1.
    # p7 (mechanical cooler): HEC_AirMech is proportional to U_MechCool -> gate by U3
    #    (as coded it is a permanent dehumidifier toward Pws(I6) = Pws(0 C)).
    # p8 (thermal screen): HEC_AirThScr = 1.7 * U_ThScr * |T_air - T_ThScr|^0.33
    #    (as coded: no U1 factor, and |T2 - VPsat_ThmScr| mixes K with Pa -> ~7x too large).
    if hasproperty(p, :vapour_fix) && p.vapour_fix > 0.5
        p7 = U3 * p7
        p8 = V1 < VPsat_ThmScr ? 0.0 : 6.4e-9 * 1.7 * U1 * (abs(T2 - T_ThmScr))^0.33 * (V1 - VPsat_ThmScr)
    end
    # --- cond_heat (default-off: identical to V2.2) ----------------------------
    # Latent heat released by screen condensation (gamma2 * p8, W m^-2) is otherwise lost.
    # The screen has no energy balance here (T_ThmScr is algebraic), so a fraction f_cond
    # of it is returned to the crop air; the rest is assumed to leave through the cover.
    if hasproperty(p, :cond_heat) && p.cond_heat > 0.5
        dT2 = dT2 + p.f_cond * gamma2 * p8 / (phi2 * rho3 * alpha5)
    end

    kappa3 = (psi1 * phi2) / (omega2 * T2)
    dV1 = (1 / kappa3) * (p1 + p2 + p3 + p4 - p5 - p6 - p7 - p8 - p9)
    # Note 79 (default-off): cover_cond = 1 -> condensation on the inside of the cover (Vanthoor MV_AirCov):
    #   MV = 6.4e-9 * 1.7 |T_air - T_cov|^0.33 * (VP_air - Psat(T_cov)) per m2 cover, if VP_air > Psat(T_cov),
    # times the cover area per floor area and the share of crop air in contact with the cover when the energy screen
    # is closed (vent_scr_perm, the same through-screen exchange share used for ventilation). The latent heat goes to
    # the cover (algebraic, lost to outside), not to the air.
    if hasproperty(p, :cover_cond) && p.cover_cond > 0.5
        VPc = Pws(I4_eff)
        if V1 > VPc
            fcov = 1.0 - U1 * (1.0 - p.vent_scr_perm)
            p10  = cov_area * fcov * 6.4e-9 * 1.7 * abs(T2 - I4_eff)^0.33 * (V1 - VPc)
            dV1 -= p10 / kappa3
        end
    end

    # --- dC1: CO2 ------------------------------------------------------------
    o1 = eta13 * U4 * lamb4 / alpha6
    o2 = U10 * psi2 / alpha6
    o3 = (U2 * phi7 / alpha6) * (C1 - I10)
    o4 = A
    o5 = (C1 - I10) * (f2 + f3 + f4 + f6CO2)

    dC1 = (1 / phi2) * (o1 + o2 + o3 - o4 - o5)

    return dT1, dT2, dC1, dV1
end
