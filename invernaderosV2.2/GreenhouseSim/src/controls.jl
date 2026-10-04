# =============================================================================
# controls.jl  --  concrete, type-stable measured-control inputs
# =============================================================================
#
# Ports 2.0's `control4inference`: read the measured control trajectories
# (U1..U12, pipe setpoint U11, pipe temperature Tpipe) from the observation
# CSV and build step-function interpolants. Returned in a CONCRETE struct so
# `c.U6(t)` is type-stable in the RHS.
#
# Knots are seconds since the first sample, so call with time measured from the
# start of the experiment (SimContext.controls_t0 handles the daily offset).

"""
    ControlInputs

Concrete container of the measured control interpolants, each callable with a
time in seconds since the experiment start:

    U1..U12   (control actuators; U11 is the pipe/air temperature setpoint)
    Tpipe     (measured heating-pipe temperature, used as input I3)
"""
struct ControlInputs{T}
    U1::T; U2::T; U3::T; U4::T; U5::T; U6::T
    U7::T; U8::T; U9::T; U10::T; U11::T; U12::T
    Tpipe::T
end

"""
    load_measured_controls(csv_path) -> ControlInputs

Read the observation CSV (columns: time, ..., U1..U12, Tpipe) and build
constant (step) interpolants indexed by seconds since the first row.
"""
function load_measured_controls(csv_path::AbstractString)
    df = DataFrame(CSV.File(csv_path))
    t = Float64.(df.time) .- Float64(df.time[1])
    CI(y) = _const_interp(Float64.(y), t)   # defined in weather.jl (included first)
    return ControlInputs(
        CI(df.U1),  CI(df.U2),  CI(df.U3),  CI(df.U4),  CI(df.U5),  CI(df.U6),
        CI(df.U7),  CI(df.U8),  CI(df.U9),  CI(df.U10), CI(df.U11), CI(df.U12),
        CI(df.Tpipe),
    )
end
