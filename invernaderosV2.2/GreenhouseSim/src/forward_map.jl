# =============================================================================
# forward_map.jl  --  SimContext, ODE RHS, and the inference forward map
# =============================================================================
#
# The forward map is the function the MCMC calls thousands of times. Speed here
# depends on: (1) no globals -- parameters flow through `p`; (2) concrete types
# -- SimContext / weather / controls are all concrete; (3) a FUNCTION BARRIER
# -- `update_params` builds the new parameter NamedTuple ONCE per proposal
# (its only dynamic step), then `solve` dispatches into a fully specialized RHS.

# --- immutable bundle of the invariant simulation inputs ---------------------
struct SimBase{P<:NamedTuple, W, C}
    params::P
    weather::W
    controls::C
    start_time::Float64      # Unix seconds at the experiment start
end

# --- per-day context passed to the ODE solver as `p` -------------------------
struct SimContext{P<:NamedTuple, W, C}
    params::P
    weather::W
    controls::C
    weather_t0::Float64      # absolute Unix time at this day's start (weather)
    controls_t0::Float64     # seconds since experiment start at this day's start
    lai::Float64             # dynamic LAI for the growth path; NaN => use weather.LAI (legacy)
end

# 5-arg constructor keeps the validated climate/inference path unchanged:
# lai defaults to NaN, and the legacy `rhs!` ignores it (uses weather.LAI).
SimContext(p, w, c, wt0, ct0) = SimContext(p, w, c, wt0, ct0, NaN)

"""
    floor_balance(T2, T5, Tdeep, p) -> (dT2_extra, dT5)

Dynamic floor/soil thermal-mass state. The floor temperature `T5` is a state
(appended LAST to each RHS's u-vector) and is passed to `climate_rhs` in place of
the old fixed soil input, so the floor's existing radiative/convective couplings
respond to a moving floor. This helper adds ONE conservative convective exchange
between air and floor (air loses exactly what the floor gains) plus a deep-soil
anchor. Parameters `h_soil, C_soil, k_deep` are calibratable; deep-soil temp = the
weather Tsoil. `climate_rhs` itself is unchanged.
"""
@inline function floor_balance(T2, T5, Tdeep, p)
    (; k_ground, C_soil, phi2, rho3, alpha5) = p
    Qaf = k_ground * (T2 - T5)                    # W m^-2: air -> floor
    kfd = (hasproperty(p, :deep_k) && p.deep_k > 0.5) ? p.k_deep : k_ground   # Note 46: own floor-deep conductance
    Qfg = kfd * (T5 - Tdeep)                      # W m^-2: floor -> deep-soil anchor
    dT2_extra = -Qaf / (phi2 * rho3 * alpha5)     # K s^-1: air loses it
    dT5 = (Qaf - Qfg) / C_soil                    # K s^-1: floor charges/discharges (fixed capacity)
    return dT2_extra, dT5
end

"""
    rhs!(du, u, ctx::SimContext, t)

Right-hand side of the greenhouse ODE system.
State `u = [T1, T2, C1, V1, integral_error, T5_floor]` (T5 = floor thermal mass).
"""
function rhs!(du, u, ctx::SimContext, t)
    @inbounds begin
        T1 = u[1]; T2 = u[2]; C1 = u[3]; V1 = u[4]; T5 = u[6]
    end
    p = ctx.params
    w = ctx.weather
    c = ctx.controls

    tcal  = t + ctx.weather_t0     # absolute time for weather interpolants
    tctrl = t + ctx.controls_t0    # experiment-relative time for controls

    # environmental inputs
    I1  = w.LAI(tcal)
    I2  = w.Iglobal(tcal)
    I4  = w.Tsky(tcal)
    I5  = w.Tout(tcal)
    I6  = w.TmechCool(tcal)
    Tdeep = w.Tsoil(tcal)
    # Note 43 (default-off): soil_dyn = 1 -> the deep-soil anchor is a model state u[9] that slowly follows the
    # floor temperature, dTdeep/dt = (T5 - Tdeep)/sd_tau, instead of the fixed weather Tsoil input.
    dTd = 0.0; dTm = 0.0
    if hasproperty(p, :soil_dyn) && p.soil_dyn > 0.5 && length(u) >= 9
        @inbounds Tdeep = u[9]
        dTd = (T5 - Tdeep) / p.sd_tau
    end
    I7  = T5                       # dynamic floor temperature (thermal mass)
    I8  = w.WindSpeed(tcal)
    (; eta1, tau1, eta2) = p
    I9  = (1 - eta1) * tau1 * eta2 * w.Idocel(tcal)
    I10 = w.CO2out(tcal)
    I11 = w.VPout(tcal)

    # photosynthesis (2.0-faithful: LAI arg = I9; see note in photosynthesis.jl)
    A = assimilation(T2, C1, I9, I9, p)

    # measured controls
    U1  = c.U1(tctrl);  U2  = c.U2(tctrl);  U3  = c.U3(tctrl);  U4 = c.U4(tctrl)
    U5  = c.U5(tctrl);  U6  = c.U6(tctrl);  U7  = c.U7(tctrl);  U8 = c.U8(tctrl)
    U9  = c.U9(tctrl);  U10 = c.U10(tctrl); U11 = c.U11(tctrl); U12 = c.U12(tctrl)
    I3  = c.Tpipe(tctrl)

    # --- Note 34 flags (default-off; need ControlInputsX and, for states, length(u) >= 8) ---
    # black_scr: blackout screen as extra night insulation, U1_eff = 1-(1-U1)(1-k_black*BlackScr)
    # vc_night_log / vc_min_night (Note 38, need vent_ctrl): at night (Iglobal < vc_night_I) replay the
    #            logged opening, and/or impose a minimum opening.
    # vent_ctrl: closed-loop roof vents. Target = clamp(vc_gain*(T2 - VentSP - vc_e0), 0, 1),
    #            raised to vc_hum*clamp((HumDefSP - HD)/vc_hband, 0, 1) when too humid; the
    #            opening is state u[7] tracking the target with time constant vc_tau [s].
    # pipe_lag:  (with heat_ctrl) the controller pipe temperature is state u[8], first-order
    #            lag pl_tau [s] toward the PI request.
    xc = hasproperty(c, :VentSP)
    if xc && hasproperty(p, :black_scr) && p.black_scr > 0.5
        U1 = 1.0 - (1.0 - U1) * (1.0 - p.k_black * c.BlackScr(tctrl))
    end
    dU8s = 0.0; dTp = 0.0
    if xc && hasproperty(p, :vent_ctrl) && p.vent_ctrl > 0.5 && length(u) >= 8
        (; vc_gain, vc_e0, vc_tau, vc_hum, vc_hband, psi1, omega2) = p
        @inbounds U8s = u[7]
        HD   = (Pws(T2) - V1) * psi1 / (omega2 * T2) * 1000.0      # model humidity deficit [g/m3]
        tgt  = clamp(vc_gain * (T2 - c.VentSP(tctrl) - vc_e0), 0.0, 1.0)
        tgt  = max(tgt, vc_hum * clamp((c.HumDefSP(tctrl) - HD) / vc_hband, 0.0, 1.0))
        dU8s = (tgt - U8s) / vc_tau
        U8   = clamp(U8s, 0.0, 1.0)
        # Note 38 night-vent options (default-off; vc_night_log = 0 and vc_min_night = 0 leave this a no-op)
        if hasproperty(p, :vc_night_log) && (p.vc_night_log > 0.5 || p.vc_min_night > 0.0) && I2 < p.vc_night_I
            if p.vc_night_log > 0.5
                U8   = clamp(c.U8(tctrl), 0.0, 1.0)      # night: replay the logged opening
                dU8s = (U8 - U8s) / vc_tau               # vent state follows the log (consistent at dawn)
            end
            U8 = max(U8, p.vc_min_night)                 # minimum night opening
        end
        # Note 42: logged minimum-vent setpoint as an input (default-off; needs ControlInputsX with VentMin)
        if hasproperty(p, :vc_min_log) && p.vc_min_log > 0.5 && hasproperty(c, :VentMin) && !(c.VentMin isa Nothing)
            U8 = max(U8, c.VentMin(tctrl))
        end
    end

    # --- heating controller (default-off: heat_ctrl = 0 -> identical to V2.2) --------
    # When the logged pipe temperature is missing/off (PipeGrow is dead in AGC-2018 AiCU),
    # a PI controller drives the pipe to hold the air at the heating setpoint U11 + hc_dT:
    #   Tpipe_ctrl = T2 + clamp(hc_Kp*e + hc_Ki*u5, 0, hc_Tmax - T2),  e = U11 + hc_dT - T2
    # The pipe used is max(input pipe, controller pipe), so it can only ADD heat. The
    # existing integral state u[5] gets anti-windup: it stays >= 0 and stops integrating
    # while the pipe is at hc_Tmax. hasproperty() is resolved at compile time for a
    # NamedTuple, so config files without these fields behave exactly as before.
    eI = U11 - T2
    if hasproperty(p, :heat_ctrl) && p.heat_ctrl > 0.5
        (; hc_Kp, hc_Ki, hc_Tmax, hc_dT) = p
        @inbounds u5 = u[5]
        e    = U11 + hc_dT - T2
        dmd  = hc_Kp * e + hc_Ki * max(u5, 0.0)        # requested pipe-air difference [K]
        room = max(hc_Tmax - T2, 0.0)
        Ttgt = T2 + clamp(dmd, 0.0, room)
        # Note 42 (default-off): hc_log_pipe = 1 -> where the logged pipe is warm (> air + 1 K) the logged pipe
        # alone heats; the PI only acts where no warm pipe is logged, and its integral is frozen otherwise.
        logwarm = hasproperty(p, :hc_log_pipe) && p.hc_log_pipe > 0.5 && I3 > T2 + 1.0
        if logwarm
            # keep the logged I3 unchanged
        elseif hasproperty(p, :pipe_lag) && p.pipe_lag > 0.5 && length(u) >= 8
            @inbounds Tps = u[8]
            dTp = (Ttgt - Tps) / p.pl_tau
            I3  = max(I3, Tps)
        else
            I3  = max(I3, Ttgt)
        end
        eI   = (logwarm || (dmd >= room && e > 0.0) || (u5 <= 0.0 && e < 0.0)) ? 0.0 : e
    end

    # Note 58 (default-off): vent_split = 1 -> leeward and windward roof vents as two flows. Inputs: U8 = leeward
    # opening (VentLee/100), U7 = windward opening (Ventwind/100; forced ventilation U7 is then 0). Each side is half
    # the roof-vent area with its own wind-pressure multiplier (vs_clee, vs_cwind) and discharge multiplier for the
    # windward side (vs_gwind); shared stack term and shared curvature vs_exp:
    #   flow = 1/2 [ lee^e sqrt(B + c_lee n3 v^2) + g_w wind^e sqrt(B + c_wind n3 v^2) ]
    # passed to climate_rhs as the single equivalent opening U8 that gives the same flow (defaults reproduce
    # U8 = (lee + wind)/2 exactly). Meant for open loop (vent_ctrl = 0) with logged openings.
    if hasproperty(p, :vent_split) && p.vent_split > 0.5
        Uw = clamp(U7, 0.0, 1.0); Ul = clamp(U8, 0.0, 1.0); U7 = 0.0
        e  = p.vs_exp
        n3 = p.nu2 * (1 - p.eta11 * U5)
        B  = p.omega1 * p.nu6 * (T2 - I5) / (T2 + I5)
        w2 = n3 * I8^2
        f0 = max(B + w2, 0.0)
        sl = sqrt(max(B + p.vs_clee * w2, 0.0)); sw = sqrt(max(B + p.vs_cwind * w2, 0.0))
        num = 0.5 * (Ul^e * sl + p.vs_gwind * Uw^e * sw)
        U8 = f0 > 1e-12 ? num / sqrt(f0) : 0.5 * (Ul^e + p.vs_gwind * Uw^e)
    end

    dT1, dT2, dC1, dV1 = climate_rhs(T1, T2, C1, V1,
        I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A,
        U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12, p)

    @inbounds begin
        dT2x, dT5 = floor_balance(T2, T5, Tdeep, p)
        # Note 46 (default-off): floor_full closes the floor energy balance; solar_split adds sun on the floor
        if (hasproperty(p, :floor_full) && p.floor_full > 0.5) || (hasproperty(p, :solar_split) && p.solar_split > 0.5)
            c1, c2, Qf = floor_full_fluxes(T1, T2, T5, I1, I2, I4, I5, U1, p)
            dT1 += c1; dT2x += c2; dT5 += Qf / p.C_soil
        end
        # Note 48 (default-off): mass_node = 1 -> lumped internal thermal mass (structure, substrate, crop water,
        # pipe water) as state u[10]: C_mass dTm/dt = h_mass (T2 - Tm) + f_mass * Iglob; the air loses h_mass (T2 - Tm).
        if hasproperty(p, :mass_node) && p.mass_node > 0.5 && length(u) >= 10
            @inbounds Tm = u[10]
            Qam  = p.h_mass * (T2 - Tm)
            dT2x -= Qam / (p.phi2 * p.rho3 * p.alpha5)
            dTm  = (Qam + p.f_mass * I2) / p.C_mass
        end
        # Note 54 (default-off): lamp_full = 1 -> the lamps release their full electrical power lamp_Pel*U12
        # (AGC-2018 ReadMe: 187 umol m-2 s-1 at 1.6 umol/J = 117 W m-2). V2.0 puts only (eta15+eta16)*alpha12 = 29 W m-2
        # into the air and (1-eta1)*alpha2*eta14*alpha12 = 4 W m-2 into the canopy (~29 % of the electrical input; the
        # lamp PAR is missing from the energy balance). Here a share lamp_fcan*(1-exp(-beta1*LAI)) heats the canopy
        # (PAR + NIR absorbed) and the rest heats the air (convection + lamp-body long-wave, lumped).
        if hasproperty(p, :lamp_full) && p.lamp_full > 0.5
            Ql   = p.lamp_Pel * U12
            # Note 58d (default-off): lamp_scr_loss > 0 -> share of the lamp power lost upward through the roof when
            # the screens are open; closure = 1 - (1 - EnScr)(1 - BlackScr). Retained = Ql * (1 - loss * (1 - closure)).
            if hasproperty(p, :lamp_scr_loss) && p.lamp_scr_loss > 0.0
                bs   = hasproperty(c, :BlackScr) && !(c.BlackScr isa Nothing) ? c.BlackScr(tctrl) : 0.0
                clos = 1.0 - (1.0 - clamp(c.U1(tctrl), 0.0, 1.0)) * (1.0 - clamp(bs, 0.0, 1.0))
                Ql  *= 1.0 - p.lamp_scr_loss * (1.0 - clos)
            end
            fc   = p.lamp_fcan * (1 - exp(-p.beta1 * I1))
            r5l  = (1 - p.eta1) * p.alpha2 * p.eta14 * p.alpha12 * U12
            h12l = (p.eta15 + p.eta16) * p.alpha12 * U12
            dT1  += (fc * Ql - r5l) / (p.alpha1 * I1)
            dT2x += ((1 - fc) * Ql - h12l) / (p.phi2 * p.rho3 * p.alpha5)
        end
        du[1] = dT1; du[2] = (dT2 + dT2x) / (hasproperty(p, :air_cap_mult) ? p.air_cap_mult : 1.0);   # Note 61: effective air heat capacity (default 1)
        du[3] = dC1; du[4] = dV1
        du[5] = eI                 # integral of the temperature error (PI); = U11 - T2 when heat_ctrl is off
        du[6] = dT5                # floor/soil thermal-mass state
    end
    nu = length(u)
    if nu > 6                      # optional states (Note 34): u[7] vent opening, u[8] pipe temperature
        @inbounds for k in 7:nu; du[k] = 0.0; end
        if nu >= 8
            @inbounds du[7] = dU8s; @inbounds du[8] = dTp
        end
        if nu >= 9
            @inbounds du[9] = dTd                # deep-soil state (Note 43); 0 when soil_dyn is off
        end
        if nu >= 10
            @inbounds du[10] = dTm               # internal thermal mass (Note 48); 0 when mass_node is off
        end
    end
    return nothing
end

# default initial state [T1, T2, C1, V1, integral_error]
const DEFAULT_U0 = [18 + 273.15, 23 + 273.15, 575.0, 1200.0, 0.0, 18 + 273.15]  # last = floor T5

"""
    simulate(base::SimBase, days, tspan, teval; u0, solver) -> Matrix

Run the day-by-day simulation with the parameters already in `base` (no
override). Returns the state matrix (rows = state, columns = time points),
carrying the end state into each following day -- the daily outer loop where a
future growth model will update LAI / biomass once per day.
"""
function simulate(base::SimBase, days, tspan, teval;
                  u0 = copy(DEFAULT_U0), solver = Tsit5())
    return _run(base.params, base, days, tspan, teval, u0, solver)
end

"""
    forward_map(x, qoi, base::SimBase, days, tspan, teval; u0, solver)
        -> (T1, T2, RH, CO2)

Override the parameters named in `qoi` with the values `x`, run the simulation,
and return the four observable trajectories used in the likelihood. This is the
function to call from the MCMC energy.
"""
function forward_map(x, qoi, base::SimBase, days, tspan, teval;
                     u0 = copy(DEFAULT_U0), solver = Tsit5())
    p = update_params(base.params, qoi, x)     # <-- the one dynamic step
    us = _run(p, base, days, tspan, teval, u0, solver)   # function barrier
    T1  = @view us[1, :]
    T2  = @view us[2, :]
    CO2 = @view us[3, :]
    RH  = rhf.(us[2, :], us[4, :])
    return T1, T2, RH, CO2
end

# Internal: type-stable daily loop. `p` is concrete, so everything below it is
# specialized by the compiler.
function _run(p::NamedTuple, base::SimBase, days, tspan, teval, u0, solver)
    us = reshape(copy(u0), :, 1)
    u_curr = copy(u0)
    for day_i in days
        ctx = SimContext(p, base.weather, base.controls,
                         base.start_time + 86400.0 * day_i,   # weather_t0
                         86400.0 * day_i)                     # controls_t0
        prob = ODEProblem(rhs!, u_curr, tspan, ctx)
        sol = solve(prob, solver; saveat = teval)
        us = hcat(us, Array(sol))
        u_curr = sol[:, end]
    end
    return us
end
