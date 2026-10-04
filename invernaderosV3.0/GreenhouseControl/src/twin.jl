# =============================================================================
# twin.jl  --  Stage 2/3: closed-loop digital twin + resource accounting
# =============================================================================
# Couples the control-climate ODE (rhs_control!) with the daily cucumber growth
# model, and now also meters the RESOURCES the controller spends each day
# (heating energy, CO2 dosed, lamp electricity) so the Stage-3 reward can price
# them. This is the environment the RL agent acts on.

"""
    day_resources(sol, setpoints, params, gT, gC) -> (heat_kWh, co2_kg, lamp_kWh_peak, lamp_kWh_off)

Recompute the controller's actuator effort over one solved day and integrate it
into resource use per m^2. Mirrors the actuator formulas in `rhs_control!`
(heating pipe power = the model's h4 term; CO2 dosing = U10*psi2/alpha6; lamps =
alpha12*U12). Peak electricity = 07:00-23:00 local.
"""
function day_resources(sol, setpoints, params, gT::PIGains, gC::PIGains)
    (; HEAT_PIPE, phi1, gamma1, n_pipes, psi2, alpha6, alpha12, eta13) = params
    FLUE_FRAC = 0.85   # fraction of boiler combustion CO2 captured+routed to the crop (free)
    ts = sol.t; nt = length(ts)
    heatW = Vector{Float64}(undef, nt); co2F = similar(heatW)
    lampW = similar(heatW); peak = falses(nt)
    @inbounds for j in 1:nt
        td = mod(ts[j], 86400.0)
        T2 = sol[2, j]; C1 = sol[3, j]; eT = sol[5, j]; eC = sol[6, j]
        Tset = setpoints.Tset(td); CO2set = setpoints.CO2_set(td); Lon = setpoints.Light_on(td)
        U_heat = clamp(gT.Kp * (Tset - T2) + gT.Ki * eT, 0.0, 1.0)
        I3  = U_heat * HEAT_PIPE + (1 - U_heat) * T2
        dTp = I3 - T2
        heatW[j] = dTp > 0 ? n_pipes * (1.99 * pi * phi1 * gamma1 * abs(dTp)^0.32) * dTp : 0.0
        U10 = clamp(gC.Kp * (CO2set - C1) + gC.Ki * eC, 0.0, 1.0)
        dose     = U10 * psi2 / alpha6             # CO2 dosing demand [mg m^-2 s^-1]
        free_co2 = FLUE_FRAC * eta13 * heatW[j]     # combustion CO2 available while heating
        co2F[j]  = max(dose - free_co2, 0.0)        # PAID (pure) CO2 only
        lampW[j] = alpha12 * Lon                  # W m^-2 electric
        h = td / 3600.0
        peak[j]  = (7.0 <= h < 23.0)
    end
    heatJ = 0.0; co2mg = 0.0; lampJp = 0.0; lampJo = 0.0
    @inbounds for j in 1:nt-1
        dt = ts[j+1] - ts[j]
        heatJ += 0.5 * (heatW[j] + heatW[j+1]) * dt
        co2mg += 0.5 * (co2F[j]  + co2F[j+1])  * dt
        lp     = 0.5 * (lampW[j] + lampW[j+1]) * dt
        (peak[j] ? (lampJp += lp) : (lampJo += lp))
    end
    return heatJ / 3.6e6, co2mg / 1e6, lampJp / 3.6e6, lampJo / 3.6e6
end

"""
    density_schedule(pairs) -> (day -> stems/m2)

Build a step schedule of stem density from (day, density) pairs, e.g.
`density_schedule([(0, 2.0), (30, 2.75), (60, 3.0)])` — planted at 2.0, interplanted
up to 2.75 at day 30 and 3.0 at day 60. Density scales the crop's per-area capacities
(veg_sink_max, set_rate, LAI_max) in `grow!`.
"""
function density_schedule(pairs::Vector{<:Tuple{<:Integer,<:Real}})
    ps = sort(collect(pairs); by = first)
    return function (d::Integer)
        v = Float64(ps[1][2])
        for (day, dens) in ps
            d >= day && (v = Float64(dens))
        end
        return v
    end
end

"""
    simulate_twin(params, gp, weather, setpoints, start_unix, ndays; ...) -> NamedTuple

Daily series: day, LAI, yield_FW [kg/m2], W_fruit, W_leaf, node, Pg, Tmean,
CO2mean, plus per-day resource use heat_kWh, co2_kg, lamp_kWh_peak, lamp_kWh_off,
and the final growth `state`.
"""
function simulate_twin(params, gp, weather, setpoints, start_unix::Float64, ndays::Int;
                       LAI0 = 0.5,
                       gT = PIGains(0.6, 0.02, 5e-4),
                       gC = PIGains(0.01, 2e-4, 5e-4),
                       solver = Tsit5(), saveat = 3600.0,
                       dtmax = 300.0,
                       density_sched = nothing,   # d -> stems/m2 (nothing = constant gp.stem_density)
                       deleaf_sched  = nothing)   # d -> LAI de-leaf target (nothing = none)
    dens0 = density_sched === nothing ? gp.stem_density : density_sched(0)
    gs = init_growth_state(gp; LAI0 = LAI0 * dens0 / 2.5)   # initial canopy scales with planting density
                                                            # (2.5 = REF_DENSITY; identity at reference)
    u  = copy(U0_CONTROL)

    day = Int[]; LAI = Float64[]; yieldFW = Float64[]; Wfruit = Float64[]; Wleaf = Float64[]
    Pg = Float64[]; Tmean = Float64[]; CO2mean = Float64[]; node = Float64[]
    heat_kWh = Float64[]; co2_kg = Float64[]; lamp_kWh_peak = Float64[]; lamp_kWh_off = Float64[]
    density = Float64[]

    for d in 0:(ndays - 1)
        u[7] = 0.0                            # reset daily assimilation integral (floor u[8] persists)
        wt0  = start_unix + 86400.0 * d
        ctx  = ControlContext(params, weather, setpoints, gs.LAI, wt0, gT, gC, params.HEAT_PIPE)
        prob = ODEProblem(rhs_control!, u, (0.0, 86400.0), ctx)
        sol  = solve(prob, solver; saveat = saveat, dtmax = dtmax)

        n    = size(sol, 2)
        Aint = sol[7, end] - sol[7, 1]
        Pgd  = Aint * 1e-3 * (30.0 / 44.0)
        Tm   = sum(@view sol[2, :]) / n - 273.15
        Cm   = sum(@view sol[3, :]) / n

        hkWh, ckg, lpk, loff = day_resources(sol, setpoints, params, gT, gC)

        _nm = String[]; _vl = Float64[]
        density_sched === nothing || (push!(_nm, "stem_density"); push!(_vl, density_sched(d)))
        deleaf_sched  === nothing || (push!(_nm, "LAI_deleaf");   push!(_vl, deleaf_sched(d)))
        gpd = isempty(_nm) ? gp : update_params(gp, _nm, _vl)
        grow!(gs, Pgd, Tm, gpd, d)
        u = sol[:, end]

        push!(day, d); push!(LAI, gs.LAI); push!(yieldFW, gs.yield_FW)
        push!(Wfruit, isempty(gs.fruits) ? 0.0 : sum(f.W for f in gs.fruits))
        push!(Wleaf, gs.W_leaf); push!(node, gs.node)
        push!(Pg, Pgd); push!(Tmean, Tm); push!(CO2mean, Cm)
        push!(heat_kWh, hkWh); push!(co2_kg, ckg)
        push!(lamp_kWh_peak, lpk); push!(lamp_kWh_off, loff)
        push!(density, density_sched === nothing ? gp.stem_density : density_sched(d))
    end

    return (; day, LAI, yield_FW = yieldFW, W_fruit = Wfruit, W_leaf = Wleaf,
             node, Pg, Tmean, CO2mean,
             heat_kWh, co2_kg, lamp_kWh_peak, lamp_kWh_off, density, state = gs)
end
