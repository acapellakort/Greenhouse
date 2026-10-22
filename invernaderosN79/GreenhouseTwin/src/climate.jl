# =============================================================================
# climate.jl -- the climate ODE (n79), one equation set, no flags
# =============================================================================
# Term names follow V2.2 (Vanthoor-style r = radiation, h = sensible heat, l = latent heat, p = vapour, o = CO2,
# f = air exchange); docs/greenhouse_model_n79.tex gives every equation with the same names.
# `model_terms` returns every term and the derivatives (used by rhs! and by the term dump of the documentation
# check); it is @inline, so rhs! keeps only what it uses and allocates nothing.

const NSTATE = 8   # [T1, T2, C1, V1, e_int, T5, Td, Tm]

struct SimContext{P<:NamedTuple, W, C}
    params::P
    weather::W
    controls::C
    weather_t0::Float64      # absolute Unix time at this day's start (weather, LAI)
    controls_t0::Float64     # seconds since the season start at this day's start (controls, crop inputs)
end

"""
    cover_temperature(T2, Tsky, Tout, Qcov, V1, U1, p) -> (Tc, fcov)

Algebraic cover node: linearised long-wave exchange kR with the crop air and the sky, convection h_cover_out to the
outside, solar Qcov, and the latent heat of the vapour condensing on it (fixed point, 4 passes).
fcov = share of the crop air in contact with the cover (1 - U1 (1 - vent_scr_perm)).
"""
@inline function cover_temperature(T2, Tsky, Tout, Qcov, V1, U1, p)
    kR  = 4.0 * p.epsil3 * p.sigma * (0.5 * (T2 + Tsky))^3
    Tc0 = (kR * (T2 + Tsky) + p.h_cover_out * Tout + Qcov) / (2.0 * kR + p.h_cover_out)
    D    = 2.0 * kR + p.h_cover_out
    fcov = 1.0 - U1 * (1.0 - p.vent_scr_perm)
    Tc   = Tc0
    for _ in 1:4
        dv = V1 - Pws(Tc)
        Tc = Tc0 + (dv > 0.0 ? p.gamma2 * fcov * 6.4e-9 * 1.7 * abs(T2 - Tc)^0.33 * dv / D : 0.0)
    end
    return Tc, fcov
end

"""
    model_terms(u, ctx, t) -> NamedTuple

All inputs, fluxes and the derivatives at state u and day-time t [s].
"""
@inline function model_terms(u, ctx::SimContext, t)
    @inbounds begin
        T1 = u[1]; T2 = u[2]; C1 = u[3]; V1 = u[4]; eint = u[5]; T5 = u[6]; Td = u[7]; Tm = u[8]
    end
    p = ctx.params; w = ctx.weather; c = ctx.controls
    tcal  = t + ctx.weather_t0
    tctrl = t + ctx.controls_t0

    # --- inputs ----------------------------------------------------------------------------------------------
    Iglob = w.Iglob(tcal); Tsky = w.Tsky(tcal); Tout = w.Tout(tcal); v = w.Wind(tcal); VPo = w.VPout(tcal)
    CO2o = w.CO2out(tcal); Tmc = w.TmechCool(tcal); L = w.LAI(tcal)
    U1 = c.EnScr(tctrl); Ubl = c.BlackScr(tctrl); U2 = c.FanPad(tctrl); U3 = c.MechCool(tctrl); U4 = c.Heater(tctrl)
    U5 = c.Shade(tctrl); U6 = c.SideVent(tctrl); U7 = c.ForcedVent(tctrl); Ulee = c.VentLee(tctrl)
    Uwin = c.VentWind(tctrl); U9 = c.Fog(tctrl); U10 = c.CO2Valve(tctrl); U11 = c.HeatSP(tctrl); U12 = c.Lamps(tctrl)
    Tlog = c.Tpipe(tctrl); Rcrop = c.CropResp(tctrl); FW = c.CropFW(tctrl)

    (; alpha1, alpha3, alpha4, alpha5, alpha6, beta1, beta3,
       delta1, delta2, delta3, delta4, delta5, delta6, delta7,
       epsil1, epsil2, epsil3, epsil4, epsil6, eps_screen,
       eta5, eta6, eta10, eta11, eta12, eta13,
       gamma, gamma1, gamma2, gamma3, gamma4, gamma5, lamb4, n_pipes,
       nu1, nu2, nu3, nu4, nu5, nu6, nu7, nu8, omega1, omega2,
       phi1, phi2, phi5, phi6, phi7, phi8, phi9, psi1, psi2, rho3, sigma, tau1, tau2, tau3,
       cover_area, vent_scr_perm, leak_cs, k_ground, k_deep, C_soil, sd_tau, C_mass, h_mass, f_mass, crop_cp,
       lamp_Pel, lamp_fcan, lamp_par_W, air_cap_mult) = p

    # --- solar column and plant light ------------------------------------------------------------------------
    so = solar_column(Iglob, U1, Ubl, L, p)
    PARp = so.ParTop + so.ParBot + lamp_par_W * clamp(U12, 0.0, 1.0)   # PAR seen by stomata and photosynthesis

    # --- photosynthesis (leaf dark respiration Q10 = 2 on the canopy temperature) -------------------------------
    rds = 2.0^((T1 - 298.15) / 10.0)
    A = assimilation(T2, C1, PARp, lai_eff(L, p.k_lai_photo), p; Rd_scale = rds)

    # --- heating: logged pipe where it is warm, PI controller to the heating setpoint elsewhere ------------------
    e    = U11 - T2
    dmd  = p.hc_Kp * e + p.hc_Ki * max(eint, 0.0)       # requested pipe-air difference [K]
    room = max(p.hc_Tmax - T2, 0.0)
    Ttgt = T2 + clamp(dmd, 0.0, room)
    logwarm = Tlog > T2 + 1.0
    Tp   = logwarm ? Tlog : max(Tlog, Ttgt)
    deint = (logwarm || (dmd >= room && e > 0.0) || (eint <= 0.0 && e < 0.0)) ? 0.0 : e   # anti-windup

    # --- roof vents: leeward and windward sides as one equivalent opening ---------------------------------------
    Ul = clamp(Ulee, 0.0, 1.0); Uw = clamp(Uwin, 0.0, 1.0)
    ve = p.vs_exp
    n3v = nu2 * (1 - eta11 * U5)
    B  = omega1 * nu6 * (T2 - Tout) / (T2 + Tout)
    w2 = n3v * v^2
    f0 = max(B + w2, 0.0)
    sl = sqrt(max(B + p.vs_clee * w2, 0.0)); sw = sqrt(max(B + p.vs_cwind * w2, 0.0))
    num = 0.5 * (Ul^ve * sl + p.vs_gwind * Uw^ve * sw)
    U8 = f0 > 1e-12 ? num / sqrt(f0) : 0.5 * (Ul^ve + p.vs_gwind * Uw^ve)

    # --- stomata and canopy transpiration conductance ------------------------------------------------------------
    q2 = Pws(T1)                                         # saturated VP at the canopy
    q7 = (1 + exp(1 / gamma5 * (PARp - delta1)))^-1      # ~1 by day, ~0 by night
    q8 = delta5 + (delta4 - delta5) * q7
    q9 = delta7 + (delta6 - delta7) * q7
    q10 = (PARp + delta2) / (PARp + delta3)
    q5 = 1 + q9 * ((q2 - V1)^2)
    q4 = 1 + q8 * ((C1 * T2 / 536.4 - 200)^2)            # CO2 in ppm
    q3 = gamma4 * q10 * q4 * q5                          # stomatal resistance
    q1 = q2 < V1 ? (2 * rho3 * alpha5 * L) / (gamma * gamma2 * gamma3) :       # dew: boundary layer only
                   (2 * rho3 * alpha5 * L) / (gamma * gamma2 * (gamma3 + q3))

    # --- air exchange rates [m3 m-2 s-1] -------------------------------------------------------------------------
    n1 = nu1 * (1 - eta10 * U5)
    n2 = nu3 * U6
    n3 = nu2 * (1 - eta11 * U5)
    f5 = n1 * n2 * v * sqrt(n3) / (2 * alpha6)          # side vents
    f6 = nu4 * sqrt(max(v, 0.25)^2 + leak_cs * abs(T2 - Tout))   # leakage (wind + stack), heat, vapour and CO2
    f7t = max(omega1 * nu6 * (T2 - Tout) / (T2 + Tout) + n3 * v^2, 0.0)
    fscr = 1.0 - U1 * (1.0 - vent_scr_perm)             # closed energy screen throttles roof ventilation
    U8e = U8 * fscr
    f7 = (U8e * nu5 * n1) / (2.0 * alpha6) * sqrt(f7t)  # roof vents
    f2 = eta6 * (U1 * f5)                                # side-vent ventilation (V2.0 form)
    f3 = U7 * phi8 / alpha6                              # forced ventilation
    f4 = eta6 * f7                                       # roof ventilation

    # --- cover node and long-wave radiation ------------------------------------------------------------------------
    Tc, fcov = cover_temperature(T2, Tsky, Tout, so.Qcov, V1, U1, p)
    g1 = 0.49 * (1 - exp(-beta3 * L))
    g2 = 1 - U1 * (1 - tau3)
    g3 = tau1 * g2 * (1 - 0.49 * pi * gamma1 * phi1) * exp(-beta3 * L)
    SkyFIR = epsil3 * sigma * (T2^4 - Tc^4)
    r6  = alpha3 * epsil1 * epsil2 * g1 * sigma * (Tp^4 - T1^4)              # pipe -> canopy
    r7  = (g1 / 0.49) * epsil2 * epsil3 * g2 * sigma * (T1^4 - Tc^4)         # canopy -> cover
    r11 = epsil4 * g3 * epsil3 * sigma * (T5^4 - Tc^4)                       # floor -> cover
    r12 = eps_screen * (U1 * tau2) * SkyFIR                                  # energy screen (air T) -> cover
    r13 = cover_area * (epsil6 * SkyFIR)                                     # air -> cover
    r14 = epsil4 * epsil2 * (1 - 0.49 * pi * gamma1 * phi1) * sigma * (T1^4 - T5^4)   # canopy -> floor
    r1  = so.Qcan                                                            # solar absorbed by the canopy

    # --- lamps: full electrical power, canopy share by absorbed light --------------------------------------------
    Ql  = lamp_Pel * U12
    fc  = lamp_fcan * (1 - exp(-beta1 * L))
    Qlc = fc * Ql
    Qla = (1 - fc) * Ql

    # --- canopy energy ---------------------------------------------------------------------------------------------
    h1 = 2 * alpha4 * L * (T1 - T2)                      # canopy -> air convection
    l1 = gamma2 * q1 * (q2 - V1)                         # transpiration latent heat
    dT1 = (r1 + Qlc + r6 - h1 - l1 - r7 - r14) / (alpha1 * L)

    # --- air energy ------------------------------------------------------------------------------------------------
    h2  = (U2 * phi7 / alpha6) * (rho3 * alpha5 * Tout - gamma2 * rho3 * (eta5 * (phi5 - phi6)))   # fan-pad
    h4  = n_pipes * (1.99 * pi * phi1 * gamma1 * abs((Tp - T2))^0.32) * (Tp - T2)            # pipe -> air
    h7  = rho3 * alpha5 * (f2 + f4 + f3 + f6) * (T2 - Tout)                                  # ventilation + leakage
    h11 = 2 * nu7 * (T2 - T5) / (phi2 + nu8)                                                 # air -> floor conduction
    l2  = gamma2 * U9 * phi9 / alpha6                                                        # fog
    Qaf = k_ground * (T2 - T5)                                                               # air -> floor convection
    Qam = h_mass * (T2 - Tm)                                                                 # air -> internal mass
    Qsa = so.Qair                                                                            # solar to the air
    Hair = h1 + h2 + h4 - h7 - h11 - r12 - r13 - l2 + Qsa - Qaf - Qam + Qla
    dT2 = Hair / (phi2 * rho3 * alpha5) / air_cap_mult

    # --- vapour ----------------------------------------------------------------------------------------------------
    T_scr = 0.5 * (Tout + T2)                            # screen temperature (algebraic)
    VPsat_scr = Pws(T_scr); VPsat_air = Pws(T2); q6 = Pws(Tmc)
    p1 = q1 * (q2 - V1)                                                                       # transpiration
    p2 = rho3 * (U2 * phi7 / alpha6) * (eta5 * (phi5 - phi6) + phi6)                          # fan-pad
    p3 = U9 * phi9 / alpha6                                                                   # fog
    p4 = eta12 * U4 * lamb4 / alpha6                                                          # direct heater
    p5 = (psi1 / omega2) * (V1 / T2 - VPo / Tout) * (f2 + f3 + f4 + f6)                       # ventilation + leakage
    p6 = (U2 * phi7 / alpha6) * (psi1 / omega2) * (V1 / T2)                                   # fan-pad outflow
    p7 = U3 * (V1 < q6 ? 0.0 : 6.4e-9 * (V1 - q6))                                            # mechanical cooling
    p8 = V1 < VPsat_scr ? 0.0 : 6.4e-9 * 1.7 * U1 * (abs(T2 - T_scr))^0.33 * (V1 - VPsat_scr)   # energy screen
    p9 = V1 < VPsat_air ? 0.0 : 6.4e-9 * 1.2 * 1000 * ((abs(T1 - T2))^0.33) * (V1 - VPsat_air) # supersaturated air
    VPc = Pws(Tc)
    p10 = V1 > VPc ? cover_area * fcov * 6.4e-9 * 1.7 * abs(T2 - Tc)^0.33 * (V1 - VPc) : 0.0   # cover
    kappa3 = (psi1 * phi2) / (omega2 * T2)
    dV1 = (p1 + p2 + p3 + p4 - p5 - p6 - p7 - p8 - p9 - p10) / kappa3

    # --- CO2 -------------------------------------------------------------------------------------------------------
    o1 = eta13 * U4 * lamb4 / alpha6                     # direct heater
    o2 = U10 * psi2 / alpha6                             # CO2 valve
    o3 = (U2 * phi7 / alpha6) * (C1 - CO2o)              # fan-pad
    o5 = (C1 - CO2o) * (f2 + f3 + f4 + f6)               # ventilation + leakage
    oR = Rcrop * 2.0^((T1 - 293.15) / 10.0)              # crop respiration (non-leaf), Q10 = 2
    dC1 = (o1 + o2 + o3 - A - o5 + oR) / phi2

    # --- floor layer, deep soil, internal mass ---------------------------------------------------------------------
    Qfd = k_deep * (T5 - Td)
    Qfs = (1.0 - f_mass) * so.Qflr
    dT5 = (Qaf - Qfd + h11 + r14 - r11 + Qfs) / C_soil
    dTd = (T5 - Td) / sd_tau
    Cm  = C_mass + crop_cp * FW
    Qms = f_mass * so.Qflr
    dTm = (Qam + Qms) / Cm

    return (; dT1, dT2, dC1, dV1, deint, dT5, dTd, dTm,
              Iglob, Tsky, Tout, v, VPo, CO2o, L, U1, Ubl, U8, U10, U11, U12, Tlog, Tp, Tc, PARp, A, so,
              q1, q3, f2, f3, f4, f6, f7, r1, r6, r7, r11, r12, r13, r14, Qlc, Qla, h1, l1, h2, h4, h7, h11, l2,
              Qaf, Qam, Qsa, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, o1, o2, o3, o5, oR, Qfd, Qfs, Qms, Cm)
end

"Right-hand side of the climate ODE, u = [T1, T2, C1, V1, e_int, T5, Td, Tm]."
function rhs!(du, u, ctx::SimContext, t)
    m = model_terms(u, ctx, t)
    @inbounds begin
        du[1] = m.dT1; du[2] = m.dT2; du[3] = m.dC1; du[4] = m.dV1
        du[5] = m.deint; du[6] = m.dT5; du[7] = m.dTd; du[8] = m.dTm
    end
    return nothing
end
