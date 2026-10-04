# =============================================================================
# growth.jl  --  cucumber source-sink dry-matter model (Marcelis-style)
# =============================================================================
#
# Daily-timestep crop growth coupled to the fast within-day climate+FvCB ODE.
# Two nested loops:
#   * within a day (seconds): climate+photosynthesis ODE at FIXED LAI, with the
#     assimilation A integrated over 24 h  (rhs_dyn! + day_assimilation);
#   * once a day: spend the day's assimilate on maintenance + organ growth, set
#     and grow fruit cohorts, harvest, prune leaves, update LAI  (grow!).
# The new LAI feeds the next day's ODE. See notes/04-cucumber-growth-model.md.
#
# Only the within-day ODE is performance-critical (and unchanged/ type-stable).
# The daily step runs ~100x per season, so it uses ordinary Julia for clarity.

# ---------------------------------------------------------------------------
# Dynamic-LAI ODE: 6 states [T1, T2, C1, V1, integral_error, integral_A]
# Uses ctx.lai (the growth model's current LAI) for BOTH the canopy (I1) and
# photosynthesis -- i.e. the real LAI, fixing the 2.0 `LAI = I9` quirk.
# ---------------------------------------------------------------------------
const K_EXT = 0.7   # canopy light-extinction coefficient (Beer's law); ~0.7 for cucumber

"Effective light-intercepting LAI (Beer's law): (1 - exp(-k*LAI))/k. Saturates as LAI grows."
@inline lai_eff(lai, k = K_EXT) = (1.0 - exp(-k * lai)) / k

# Reference planting density [stems/m2] at which the per-m2 growth parameters
# (veg_sink_max, set_rate, LAI_max) are defined. Actual `gp.stem_density` scales
# those per-area capacities by (stem_density / REF_DENSITY); at the default 2.5
# the factor is 1, so existing calibrations reproduce exactly (backward compatible).
const REF_DENSITY = 2.5

function rhs_dyn!(du, u, ctx::SimContext, t)
    @inbounds begin
        T1 = u[1]; T2 = u[2]; C1 = u[3]; V1 = u[4]; T5 = u[7]
    end
    p = ctx.params
    w = ctx.weather
    c = ctx.controls
    LAI = ctx.lai

    tcal  = t + ctx.weather_t0
    tctrl = t + ctx.controls_t0

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
    I1  = LAI                          # real leaf area index

    A = assimilation(T2, C1, I9, lai_eff(LAI), p)   # canopy light interception (Beer's law)

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
        du[5] = U11 - T2
        du[6] = A                      # integral of assimilation over the day
        du[7] = dT5                    # floor/soil thermal-mass state
    end
    return nothing
end

"""
    day_assimilation(p, weather, controls, lai, weather_t0, controls_t0, u0, tspan, teval; solver)
        -> (Pg_CH2O, Tmean_C, u_end)

Run one day of the climate+FvCB ODE at fixed `lai` and return the daily gross
assimilate `Pg` [g CH2O m^-2 d^-1], the mean air temperature [deg C], and the
end state (to carry thermal continuity into the next day). `u0` is 6-long with
u0[6]=0 (the A-integral accumulator).
"""
function day_assimilation(p::NamedTuple, weather, controls, lai,
                          weather_t0, controls_t0, u0, tspan, teval; solver = Tsit5())
    ctx = SimContext(p, weather, controls, weather_t0, controls_t0, lai)
    prob = ODEProblem(rhs_dyn!, u0, tspan, ctx)
    sol = solve(prob, solver; saveat = teval)
    Aint  = sol[6, end] - sol[6, 1]                 # mg CO2 m^-2 over the day
    Tmean = sum(@view sol[2, :]) / size(sol, 2) - 273.15   # deg C
    Pg    = Aint * 1e-3 * (30.0 / 44.0)             # -> g CH2O m^-2 d^-1
    return Pg, Tmean, sol[:, end]
end

# ---------------------------------------------------------------------------
# Growth-model parameters (cucumber literature defaults; all calibratable).
# ---------------------------------------------------------------------------
"""
    load_growth_params(json_path) -> NamedTuple

Read a flat `{name: value}` (or `{name: {val: value}}`) JSON into a concrete
Float64 NamedTuple, so growth parameters can be inferred with `update_params`
just like the climate/crop constants.
"""
function load_growth_params(json_path::AbstractString)
    raw = JSON.parsefile(json_path)
    d = Dict{Symbol,Float64}()
    for (k, v) in raw
        d[Symbol(k)] = v isa Dict ? Float64(v["val"]) : Float64(v)
    end
    return (; d...)
end

# ---------------------------------------------------------------------------
# Crop state (daily). Cohort-based fruits; ordinary mutable Julia (not hot).
# ---------------------------------------------------------------------------
mutable struct Fruit
    dev::Float64    # development stage 0..1
    W::Float64      # cohort dry weight [g/m2]
    n::Float64      # number of fruits in cohort [fruits/m2]
end

mutable struct GrowthState
    W_leaf::Float64
    W_stem::Float64
    W_root::Float64
    LAI::Float64
    fruits::Vector{Fruit}
    yield_DW::Float64      # cumulative harvested fruit dry weight [g/m2]
    yield_FW::Float64      # cumulative harvested fruit fresh weight [kg/m2]
    node::Float64          # main-stem node number
end

"""
    init_growth_state(gp; LAI0, stem0, root0) -> GrowthState

Young transplant: leaf mass set from LAI0 via SLA, small stem/root, no fruit.
"""
function init_growth_state(gp; LAI0 = 0.5, stem0 = 8.0, root0 = 4.0)
    W_leaf0 = LAI0 / gp.SLA
    return GrowthState(W_leaf0, stem0, root0, LAI0, Fruit[], 0.0, 0.0, 0.0)
end

# smoothstep S(x)=x^2(3-2x); its derivative 6x(1-x) integrates to 1 over [0,1],
# so per-fruit potential growth `Wf_max * 6 dev (1-dev) * ddev` totals Wf_max.
@inline _sink_shape(dev) = 6.0 * dev * (1.0 - dev)

"""
    grow!(s::GrowthState, Pg, Tmean_C, gp, day)

Advance the crop state one day given daily gross assimilate `Pg` [g CH2O m^-2],
mean temperature `Tmean_C`, and the season `day` index (for fruit-set timing).
Returns the total fruit dry weight on the plant.
"""
function grow!(s::GrowthState, Pg::Float64, Tmean_C::Float64, gp, day::Int)
    dT = max(Tmean_C - gp.T_base, 0.0)          # thermal drive above base

    # 1. maintenance respiration (g CH2O/m2/day)
    W_fruit_tot = isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits)
    q10f = gp.Q10_resp ^ ((Tmean_C - 25.0) / 10.0)
    Rm = q10f * (gp.k_m_leaf * s.W_leaf + gp.k_m_stem * s.W_stem +
                 gp.k_m_root * s.W_root + gp.k_m_fruit * W_fruit_tot)

    # 2. assimilate available for growth -> potential dry-matter increment
    Cav = max(Pg - Rm, 0.0)
    dDM_pot = Cav / gp.asrq                       # g DM/m2/day the source can supply

    # 3. development: node appearance stops at topping; fruit cohorts always age
    active = day < gp.top_day                       # crop head still growing?
    if active
        s.node += gp.node_rate * dT
    end
    ddev = dT / gp.DD_fruit                          # dev increment per day
    for f in s.fruits
        f.dev += ddev
    end

    # 4. sink strengths of EXISTING cohorts (potential growth, g DM/m2/day).
    # Per-area capacities scale with planting density (per-m2 = per-stem-basis x dens).
    dens    = gp.stem_density / REF_DENSITY
    lai_cap = gp.LAI_max * dens
    veg_sink = max(gp.veg_sink_max * dens * (1.0 - s.LAI / lai_cap), 0.0)
    fruit_sink_pre = isempty(s.fruits) ? 0.0 :
                     sum(f.n * gp.Wf_max * _sink_shape(f.dev) * ddev for f in s.fruits)
    total_sink_pre = veg_sink + fruit_sink_pre

    # 5. fruit set REGULATED by carbohydrate status (supply/demand feedback).
    # When the standing fruit load is carbon-limiting (supply << demand) set is
    # suppressed; after those fruits are harvested the load drops and set resumes.
    # With the ~2-week fruit-growth delay this produces the harvest "flushes".
    if active && day >= gp.set_start_day
        ss = dDM_pot / (total_sink_pre + 1e-9)       # supply / demand
        set_frac = clamp((ss - gp.set_r_low) / (gp.set_r_high - gp.set_r_low), 0.0, 1.0)
        push!(s.fruits, Fruit(0.0, 0.0, gp.set_rate * dens * set_frac))
    end

    # 6. allocate: source- OR sink-limited, split by sink share (new cohort dev=0 => 0 sink)
    fruit_sinks = [f.n * gp.Wf_max * _sink_shape(f.dev) * ddev for f in s.fruits]
    total_sink  = veg_sink + (isempty(fruit_sinks) ? 0.0 : sum(fruit_sinks))
    if total_sink > 0.0
        dDM = min(dDM_pot, total_sink)               # sink-limited if source is ample
        veg_DM = dDM * veg_sink / total_sink
        s.W_leaf += veg_DM * gp.frac_leaf
        s.W_stem += veg_DM * gp.frac_stem
        s.W_root += veg_DM * gp.frac_root
        for (i, f) in enumerate(s.fruits)
            f.W += dDM * fruit_sinks[i] / total_sink
        end
    end

    # 7. harvest mature fruit
    keep = Fruit[]
    for f in s.fruits
        if f.dev >= 1.0
            s.yield_DW += f.W
            s.yield_FW += (f.W / gp.DMC_fruit) / 1000.0   # g DW -> kg FW / m2
        else
            push!(keep, f)
        end
    end
    s.fruits = keep

    # 8. leaf pruning: physical density cap AND the management DE-LEAF target
    # (LAI_deleaf), whichever is lower. De-leafing sheds the lower, shaded leaves
    # -- which fix almost nothing (lai_eff is saturated) but still respire -- WITHOUT
    # touching stem density or fruit set. Default gp.LAI_deleaf is large (no de-leafing).
    cap = min(lai_cap, gp.LAI_deleaf)
    if gp.SLA * s.W_leaf > cap
        s.W_leaf = cap / gp.SLA
    end
    s.LAI = gp.SLA * s.W_leaf

    return W_fruit_tot
end

# ---------------------------------------------------------------------------
# Season driver: couples the two loops over `ndays`.
# ---------------------------------------------------------------------------
"""
    simulate_growth(p, gp, weather, controls, start_unix, ndays;
                    LAI0, tspan, teval, u0) -> NamedTuple of daily time series

Runs the full digital twin: each day integrates the climate+FvCB ODE at the
current LAI to get gross assimilate, then advances the cucumber growth model.
Returns daily vectors: day, LAI, W_leaf, W_stem, W_fruit, yield_FW, Pg, Tmean.
"""
function simulate_growth(p::NamedTuple, gp, weather, controls, start_unix::Float64, ndays::Int;
                         LAI0 = 0.5,
                         tspan = (0.0, 86400.0),
                         teval = 0.0:600.0:86400.0,
                         u0 = [18 + 273.15, 23 + 273.15, 575.0, 1200.0, 0.0, 0.0, 18 + 273.15])
    s = init_growth_state(gp; LAI0 = LAI0)
    u = copy(u0)

    day      = Int[]
    LAI      = Float64[]
    W_leaf   = Float64[]
    W_stem   = Float64[]
    W_fruit  = Float64[]
    yield_FW = Float64[]
    Pg_hist  = Float64[]
    Tmean_h  = Float64[]

    for d in 0:(ndays - 1)
        u[6] = 0.0                              # reset the daily A-integral
        wt0 = start_unix + 86400.0 * d
        ct0 = 86400.0 * d
        Pg, Tmean, u = day_assimilation(p, weather, controls, s.LAI, wt0, ct0, u, tspan, teval)
        grow!(s, Pg, Tmean, gp, d)

        push!(day, d)
        push!(LAI, s.LAI)
        push!(W_leaf, s.W_leaf)
        push!(W_stem, s.W_stem)
        push!(W_fruit, isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits))
        push!(yield_FW, s.yield_FW)
        push!(Pg_hist, Pg)
        push!(Tmean_h, Tmean)
    end

    return (; day, LAI, W_leaf, W_stem, W_fruit, yield_FW, Pg = Pg_hist, Tmean = Tmean_h, state = s)
end
