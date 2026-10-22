# =============================================================================
# simulate.jl -- day-by-day integration
# =============================================================================

"""
    simdays(p, season, days, u0, reltol; tspan, teval) -> Matrix or nothing

Integrate the climate ODE day by day (Tsit5, explicit relative tolerance, default absolute tolerance), carrying the
end state of each day into the next. `season` provides `weather`, `controls` and `t0` (absolute Unix time of day 0).
Returns the states as columns: the initial state, then the `teval` points of every day (nothing if a day fails).
"""
function simdays(p::NamedTuple, season, days, u0, reltol;
                 tspan = (0.0, 86400.0), teval = 0.0:300.0:86400.0)
    npd = length(teval)
    us = Matrix{Float64}(undef, length(u0), 1 + length(days) * npd); us[:, 1] = u0
    u = copy(u0); col = 2
    for d in days
        ctx = SimContext(p, season.weather, season.controls, season.t0 + 86400.0 * d, 86400.0 * d)
        ts  = collect(0.0:300.0:86400.0)   # stop at every 5-min input-step boundary
        sol = try
            solve(ODEProblem(rhs!, u, tspan, ctx), Tsit5(); saveat = collect(teval),
                  tstops = ts, reltol = reltol)
        catch
            return nothing
        end
        length(sol.t) == npd || return nothing
        S = Array(sol); us[:, col:col+npd-1] = S; col += npd; u = S[:, end]
    end
    return us
end
