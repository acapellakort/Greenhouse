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
    Qfg = k_ground * (T5 - Tdeep)                 # W m^-2: floor -> deep-soil anchor
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

    dT1, dT2, dC1, dV1 = climate_rhs(T1, T2, C1, V1,
        I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A,
        U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12, p)

    @inbounds begin
        dT2x, dT5 = floor_balance(T2, T5, Tdeep, p)
        du[1] = dT1; du[2] = dT2 + dT2x; du[3] = dC1; du[4] = dV1
        du[5] = U11 - T2           # integral of the temperature error (PI)
        du[6] = dT5                # floor/soil thermal-mass state
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
