# =============================================================================
# inputs.jl -- concrete, type-stable model inputs
# =============================================================================
# Every input is a callable of time. Weather inputs and the LAI are called with ABSOLUTE time (Unix seconds);
# controls, crop respiration and crop fresh mass with time since the season start (SimContext gives both).

"Step (zero-order hold) function on knots t with constant extrapolation: v[k] for t[k] <= x < t[k+1]."
struct StepFun
    t::Vector{Float64}
    v::Vector{Float64}
end
@inline (f::StepFun)(x) = @inbounds f.v[clamp(searchsortedlast(f.t, x), 1, length(f.v))]

"A constant input."
struct ConstFun
    v::Float64
end
@inline (c::ConstFun)(_x) = c.v

"""
    WeatherInputs

Outside conditions and the crop LAI (absolute time):
    Iglob [W m-2] global radiation, Tsky [K] sky temperature, Tout [K], Wind [m s-1], VPout [Pa],
    CO2out [mg m-3], TmechCool [K] mechanical-cooling surface temperature, LAI [m2 leaf m-2 floor].
"""
struct WeatherInputs{A,B,C,D,E,F,G,H}
    Iglob::A
    Tsky::B
    Tout::C
    Wind::D
    VPout::E
    CO2out::F
    TmechCool::G
    LAI::H
end

"""
    ControlInputs

Actuators and crop inputs (time since season start). Fractions 0-1 unless stated.
    EnScr (U1) energy screen | BlackScr blackout screen | FanPad (U2) | MechCool (U3) | Heater (U4) direct air heater
    Shade (U5) | SideVent (U6) | ForcedVent (U7) | VentLee, VentWind roof vents (leeward, windward) | Fog (U9)
    CO2Valve (U10) | HeatSP (U11) [K] heating setpoint | Lamps (U12) | Tpipe [K] logged pipe temperature
    CropResp [mg CO2 m-2 s-1] total crop respiration at 20 C (kept for backward compat) | CropFW [kg m-2] standing
    fresh mass of stems, roots and fruits on the plant.
    CropRm20 [mg CO2 m-2 s-1] maintenance respiration at 20 C (Q10 applied live to T1) |
    CropRg [mg CO2 m-2 s-1] growth respiration (no extra Q10, driven by carbon availability).
"""
struct ControlInputs{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19}
    EnScr::T1
    BlackScr::T2
    FanPad::T3
    MechCool::T4
    Heater::T5
    Shade::T6
    SideVent::T7
    ForcedVent::T8
    VentLee::T9
    VentWind::T10
    Fog::T11
    CO2Valve::T12
    HeatSP::T13
    Lamps::T14
    Tpipe::T15
    CropResp::T16
    CropFW::T17
    CropRm20::T18   # maintenance respiration at 20 C [mg CO2 m-2 s-1], Q10 applied live
    CropRg::T19     # growth respiration [mg CO2 m-2 s-1], no extra Q10
end
