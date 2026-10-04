# =============================================================================
# weather.jl  --  concrete, type-stable weather inputs
# =============================================================================
#
# Ports 2.0's `readweatherdata_NoUnits.jl`, but returns a CONCRETE struct of
# callable interpolants instead of loose globals / a Dict{Symbol,Any}. Every
# field is called in the RHS as `w.Tout(t)`; because the struct is concrete,
# each such call is type-stable.

"A zero-cost callable that ignores its argument and returns a constant."
struct Const{T}
    v::T
end
(c::Const)(_t) = c.v

# Version-agnostic constant (step) interpolation with constant extrapolation.
# DataInterpolations changed its API: older versions take `extrapolate=true`,
# newer ones take `extrapolation=ExtrapolationType.Constant`. Try new, fall back.
function _const_interp(u, t)
    try
        return ConstantInterpolation(u, t; extrapolation = DataInterpolations.ExtrapolationType.Constant)
    catch
        return ConstantInterpolation(u, t; extrapolate = true)
    end
end

"""
    WeatherInputs

Concrete container of the environmental driving functions. Fields are callable
with an absolute time (Unix seconds):

    Iglobal, Tsky, Tout, TmechCool, Tsoil, WindSpeed, CO2out, VPout, Idocel, LAI
"""
struct WeatherInputs{TI,TK,TT,TM,TS,TW,TC,TV,TD,TL}
    Iglobal::TI      # global / external shortwave radiation      [W m^-2]
    Tsky::TK         # sky temperature                            [K]
    Tout::TT         # outside air temperature                    [K]
    TmechCool::TM    # mechanical-cooling temperature             [K]
    Tsoil::TS        # soil temperature                           [K]
    WindSpeed::TW    # wind speed                                 [m s^-1]
    CO2out::TC       # outside CO2 concentration                  [mg m^-3]
    VPout::TV        # outside vapour pressure                    [Pa]
    Idocel::TD       # radiation reaching the cover (= Iglobal)   [W m^-2]
    LAI::TL          # leaf area index (input; growth model TBD)  [-]
end

# --- saturation-vapour-pressure helper used to derive outside VP from RH -----
_q2(T) = T > 0 ?
    611.21 * exp((18.678 - (T / 234.5)) * (T / (257.14 + T))) :
    611.21 * exp((23.036 - (T / 333.7)) * (T / (279.82 + T)))
_vapour_pressure(T_celsius, RH) = _q2(T_celsius) * (RH / 100.0)

"""
    _sky_temp(Tc, RH, cloud_pct) -> Tsky [K]

Effective sky (radiative) temperature from air temperature [C], relative humidity
[%] and total cloud cover [%]. Berdahl-Martin clear-sky emissivity from dew point,
raised toward 1 (Tsky -> Tout) as cloud cover increases. Replaces the unphysical
fixed cold sky that made the twin over-radiate at night.
"""
@inline function _sky_temp(Tc::Float64, RH::Float64, cloud_pct::Float64)
    isnan(Tc) && return -0.4 + 273.15
    rh = clamp(isnan(RH) ? 80.0 : RH, 1.0, 100.0)
    a = 17.27; b = 237.7
    g = (a*Tc)/(b+Tc) + log(rh/100)
    Tdew = b*g/(a-g)
    eps_clear = clamp(0.711 + 0.56*(Tdew/100) + 0.73*(Tdew/100)^2, 0.6, 1.0)
    c = clamp((isnan(cloud_pct) ? 0.0 : cloud_pct)/100, 0.0, 1.0)
    eps_sky = eps_clear + (1 - eps_clear)*c
    return (Tc + 273.15) * eps_sky^0.25
end

"""
    load_weather(csv_path, start_date, end_date;
                 Tsky, TmechCool, Tsoil, CO2out, LAI) -> WeatherInputs

Load and filter the meteo CSV (same layout as 2.0's `dataset_meteo_holanda.csv`:
col 1 = time, 7 = temperature °C, 8 = RH %, 13 = shortwave radiation, 14 = wind
km/h) between `start_date` and `end_date`, and build interpolants indexed by
absolute Unix time. Sky / mech-cool / soil temperatures, outside CO2 and LAI
default to the constant values used in 2.0 but are keyword-overridable.
"""
function prepare_meteo(csv_path::AbstractString)
    df = CSV.read(csv_path, DataFrame)
    cloudcol = [ (x isa Number ? Float64(x) : NaN) for x in df[!, 9] ]
    rename!(df, 1  => :Unix_Time, 7  => :Temperature, 8 => :Relative_Humidity,
                9  => :Vapor_Pressure, 13 => :External_Radiation, 14 => :Wind_Speed)
    df.Unix_Time = [datetime2unix(DateTime(r.Year, r.Month, r.Day, r.Hour, r.Minute))
                    for r in eachrow(df)]
    df.Vapor_Pressure = _vapour_pressure.(df.Temperature, df.Relative_Humidity)
    df[!, :CloudPct]  = cloudcol
    return df
end

"Build WeatherInputs for a [start,end] window from an already-prepared meteo df."
function weather_from_df(df, start_date::DateTime, end_date::DateTime;
                         Tsky::Float64 = -0.4 + 273.15, TmechCool::Float64 = 273.15,
                         Tsoil::Float64 = 18.0 + 273.15, CO2out::Float64 = 834.7,
                         LAI::Float64 = 2.0, sky_from_clouds::Bool = false)
    u0 = datetime2unix(start_date); u1 = datetime2unix(end_date)
    sub = df[(df.Unix_Time .>= u0) .& (df.Unix_Time .<= u1), :]
    nrow(sub) > 1 || error("weather_from_df: window $start_date .. $end_date has $(nrow(sub)) rows " *
                           "(gap in the meteo record?). Pick a covered year.")
    t     = Float64.(sub.Unix_Time)
    Tout  = Float64.(sub.Temperature) .+ 273.15
    VPout = Float64.(sub.Vapor_Pressure)
    Iglob = Float64.(sub.External_Radiation)
    Wind  = Float64.(sub.Wind_Speed) .* (5.0 / 18.0)
    CI(y) = _const_interp(y, t)
    Ig = CI(Iglob)
    Tsky_field = if sky_from_clouds
        Tc = Float64.(sub.Temperature); RHs = Float64.(sub.Relative_Humidity); CC = Float64.(sub.CloudPct)
        CI(_sky_temp.(Tc, RHs, CC))
    else
        Const(Tsky)
    end
    return WeatherInputs(Ig, Tsky_field, CI(Tout), Const(TmechCool), Const(Tsoil),
                         CI(Wind), Const(CO2out), CI(VPout), Ig, Const(LAI))
end

function load_weather(csv_path::AbstractString, start_date::DateTime, end_date::DateTime; kwargs...)
    return weather_from_df(prepare_meteo(csv_path), start_date, end_date; kwargs...)
end