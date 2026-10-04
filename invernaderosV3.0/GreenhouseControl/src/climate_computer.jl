# =============================================================================
# climate_computer.jl  --  setpoint + PID/proportional control layer (Stage 1)
# =============================================================================
#
# Mimics a greenhouse environmental computer: the agent (or a grower) programs
# SETPOINT schedules, and low-level controllers translate them into the actuators
# U1..U12 that drive the V2.2 climate ODE. Control is on greenhouse AIR temp T2
# (the physically-measured, controlled variable), not canopy T1.
#
# Controllers (ported/cleaned from GreenHouse_3.0):
#   - heating:      continuous leaky-integral PI on T2 -> pipe temperature I3
#   - CO2 dosing:   continuous leaky-integral PI on C1 -> valve U10
#   - ventilation:  proportional band on T2 -> side U6 / roof U8
#   - thermal screen U1: rule (outside temp < ToutMax), lights U12: schedule
# PI integrals are ODE states (u[5]=temp, u[6]=CO2), so the controller has memory
# consistent with the solver.

# --- zero-cost constant callable ---------------------------------------------
struct _Const{T}; v::T; end
(c::_Const)(_t) = c.v

# version-agnostic constant (step) interpolation with constant extrapolation
function _const_interp(u, t)
    try
        return ConstantInterpolation(u, t; extrapolation = DataInterpolations.ExtrapolationType.Constant)
    catch
        return ConstantInterpolation(u, t; extrapolate = true)
    end
end

"Build a within-day step schedule from hour marks and values (padded to [0,24h])."
function _sched(hours::Vector{<:Real}, vals::Vector{<:Real})
    t = Float64.(hours) .* 3600.0; v = Float64.(vals)
    if t[1]   > 0.0;     t = vcat(0.0, t);      v = vcat(v[1], v);   end
    if t[end] < 86400.0; t = vcat(t, 86400.0);  v = vcat(v, v[end]); end
    return _const_interp(v, t)
end

# --- setpoint schedules (the agent's action surface) -------------------------
"Setpoint schedules, each callable with within-day seconds (0..86400)."
struct Setpoints{A,B,C,D,E,F,G,H}
    Tset::A         # air-temperature setpoint [K]
    CO2_set::B      # CO2 setpoint [mg/m3]
    Light_on::C     # lamp on/off (0/1)
    VentpBand::D    # ventilation proportional band [K]
    ofset::E        # ventilation dead-band above setpoint [K]
    ToutMax::F      # outside-temp threshold: also close the energy screen when colder [K]
    ScreenRad::G    # radiation below which the energy screen closes at night [W/m2]
    VPDmin::H       # vapour-pressure-deficit floor [kPa]: below it (too humid), dehumidify (vent + crack screen)
end

"""
    daynight_setpoints(; ...) -> Setpoints

Convenience builder for a day/night climate program. The RL agent will later set
these values (and/or richer hourly schedules) directly.
"""
function daynight_setpoints(; Tset_day = 294.15, Tset_night = 291.15,
                             day_start = 6.0, day_end = 20.0,
                             CO2_set = 800.0,
                             light_start = 6.0, light_end = 20.0,
                             VentpBand = 5.0, ofset = 1.0, ToutMax = 295.0,
                             ScreenRad = 50.0, VPDmin = 0.5, lamp_intensity = 1.0)
    Tsched = _sched([0.0, day_start, day_end],   [Tset_night, Tset_day, Tset_night])
    Lsched = _sched([0.0, light_start, light_end], [0.0, lamp_intensity, 0.0])   # 0..1 dimmable
    return Setpoints(Tsched, _Const(CO2_set), Lsched,
                     _Const(VentpBand), _Const(ofset), _Const(ToutMax),
                     _Const(ScreenRad), _Const(VPDmin))
end

# --- controllers -------------------------------------------------------------
"PI gains for a continuous leaky-integral controller."
struct PIGains
    Kp::Float64; Ki::Float64; decay::Float64
end

"Proportional-band ventilation on air temp -> (U_side, U_roof), each in [0,1]."
function vent_control(Tset, T2, pBand, ofset)
    if T2 < Tset + ofset
        Uvent = 0.0
    elseif T2 < Tset + pBand + ofset
        Uvent = (T2 - Tset - ofset) / pBand
    else
        Uvent = 1.0
    end
    return 2 * min(Uvent, 0.5), 2 * max(0.0, Uvent - 0.5)   # side, roof
end

"""
    screen_control(I2, I5, VPD, ScreenRad, ToutMax, VPDmin) -> U1 in [0,1]

Energy screen: closed (1) at night (radiation I2 < ScreenRad) or when outside is
cold (I5 < ToutMax) to retain heat; cracked to 0.85 when it gets too humid
(VPD < VPDmin) to release moisture; open (0) by day. VPD in kPa.
"""
@inline function screen_control(I2, I5, VPD, ScreenRad, ToutMax, VPDmin)
    closed = (I2 < ScreenRad || I5 < ToutMax) ? 1.0 : 0.0
    return (closed > 0.0 && VPD < VPDmin) ? 0.85 : closed
end

"Extra ventilation fraction [0,1] for dehumidification when VPD falls below VPDmin (too humid). VPD, VPDmin in kPa."
@inline function humidity_vent(VPD, VPDmin; band = 0.2, vmax = 0.6)
    VPD >= VPDmin ? 0.0 : clamp((VPDmin - VPD)/band, 0.0, 1.0) * vmax
end

# --- simulation context passed to the ODE as `p` -----------------------------
struct ControlContext{P<:NamedTuple, W, S}
    params::P
    weather::W
    setpoints::S
    lai::Float64          # fixed within this run (growth coupling = Stage 2)
    weather_t0::Float64   # absolute Unix time at t=0
    gT::PIGains
    gC::PIGains
    heat_pipe::Float64    # max heating-pipe temperature [K]
end

"""
    rhs_control!(du, u, ctx::ControlContext, t)

Closed-loop RHS. State u = [T1, T2, C1, V1, tempInteg, co2Integ].
"""
function rhs_control!(du, u, ctx::ControlContext, t)
    @inbounds begin
        T1 = u[1]; T2 = u[2]; C1 = u[3]; V1 = u[4]; eT = u[5]; eC = u[6]; T5 = u[8]
    end
    p = ctx.params; w = ctx.weather; sp = ctx.setpoints
    tcal = t + ctx.weather_t0
    td   = mod(t, 86400.0)                     # within-day time for setpoints

    I2 = w.Iglobal(tcal); I4 = w.Tsky(tcal); I5 = w.Tout(tcal); I6 = w.TmechCool(tcal)
    Tdeep = w.Tsoil(tcal); I7 = T5; I8 = w.WindSpeed(tcal)  # I7 = dynamic floor (thermal mass)
    (; eta1, tau1, eta2) = p
    I9 = (1 - eta1) * tau1 * eta2 * w.Idocel(tcal)
    I10 = w.CO2out(tcal); I11 = w.VPout(tcal); I1 = ctx.lai

    # setpoints
    Tset = sp.Tset(td); CO2set = sp.CO2_set(td); Lon = sp.Light_on(td)
    pBand = sp.VentpBand(td); ofs = sp.ofset(td); ToutMax = sp.ToutMax(td)
    ScreenRad = sp.ScreenRad(td); VPDmin = sp.VPDmin(td)
    VPD = (Pws(T2) - V1) / 1000.0            # vapour-pressure deficit [kPa]

    # controllers (on air temp T2 and CO2 C1)
    U1 = screen_control(I2, I5, VPD, ScreenRad, ToutMax, VPDmin)   # energy screen
    U6, U8 = vent_control(Tset, T2, pBand, ofs)
    Uhv = humidity_vent(VPD, VPDmin); U6 = max(U6, Uhv); U8 = max(U8, Uhv)  # dehumidify
    U12 = Lon
    errT   = Tset - T2
    U_heat = clamp(ctx.gT.Kp * errT + ctx.gT.Ki * eT, 0.0, 1.0)
    I3     = U_heat * ctx.heat_pipe + (1 - U_heat) * T2
    errC = CO2set - C1
    U10  = clamp(ctx.gC.Kp * errC + ctx.gC.Ki * eC, 0.0, 1.0)
    U2 = U3 = U4 = U5 = U7 = U9 = 0.0

    A = assimilation(T2, C1, I9, GreenhouseSim.lai_eff(I1), p; VPD = VPD)

    dT1, dT2, dC1, dV1 = climate_rhs(T1, T2, C1, V1,
        I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A,
        U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12, p)

    # anti-windup: freeze the integral when the actuator is saturated and the error
    # would push it further into saturation (prevents the CO2/heating runaway that
    # windup causes once the leakage is realistic).
    satT = (U_heat >= 1.0 && errT > 0.0) || (U_heat <= 0.0 && errT < 0.0)
    satC = (U10    >= 1.0 && errC > 0.0) || (U10    <= 0.0 && errC < 0.0)

    @inbounds begin
        dT2x, dT5 = floor_balance(T2, T5, Tdeep, p)
        du[1] = dT1; du[2] = dT2 + dT2x; du[3] = dC1; du[4] = dV1
        du[5] = satT ? (-ctx.gT.decay * eT) : (errT - ctx.gT.decay * eT)   # PI + anti-windup (temp)
        du[6] = satC ? (-ctx.gC.decay * eC) : (errC - ctx.gC.decay * eC)   # PI + anti-windup (CO2)
        du[7] = A                              # integral of assimilation over the day
        du[8] = dT5                            # floor/soil thermal-mass state
    end
    return nothing
end

# default initial state [T1, T2, C1, V1, tempInteg, co2Integ, assimInteg]
const U0_CONTROL = [18 + 273.15, 21 + 273.15, 800.0, 1200.0, 0.0, 0.0, 0.0, 18 + 273.15]  # last = floor T5

"""
    run_control(params, weather, setpoints, start_unix, ndays; lai, u0, gT, gC, saveat)

Simulate `ndays` under a setpoint program at fixed LAI. Returns the ODE solution
(state u = [T1,T2,C1,V1,tempInteg,co2Integ]).
"""
function run_control(params, weather, setpoints, start_unix::Float64, ndays::Int;
                     lai = 2.0, u0 = copy(U0_CONTROL),
                     gT = PIGains(0.2, 0.01, 0.001),
                     gC = PIGains(0.005, 1e-4, 0.001),
                     solver = Tsit5(), saveat = 600.0)
    ctx  = ControlContext(params, weather, setpoints, lai, start_unix, gT, gC, params.HEAT_PIPE)
    prob = ODEProblem(rhs_control!, u0, (0.0, ndays * 86400.0), ctx)
    return solve(prob, solver; saveat = saveat)
end
