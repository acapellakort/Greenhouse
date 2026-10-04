module Controls

using DataInterpolations
using JSON
using Dates
using DataFrames
using CSV

export ControlSet, MeasuredControls, load_controls, load_control_interpolations, get_val

# --- 1. JSON Date-Aware Controls ---
struct ControlSet
    data::Dict{Symbol, Vector{Pair{Date, ConstantInterpolation}}}
end

function load_controls(json_path::String)
    raw = JSON.parsefile(json_path)
    data = Dict{Symbol, Vector{Pair{Date, ConstantInterpolation}}}()
    
    for (var_name, dates_dict) in raw
        var_sym = Symbol(var_name)
        date_list = Pair{Date, ConstantInterpolation}[]
        for (date_str, time_dict) in dates_dict
            d = Date(date_str, "dd-mm-yyyy")
            push!(date_list, d => build_daily_interpolation(time_dict))
        end
        sort!(date_list, by=x->x.first)
        data[var_sym] = date_list
    end
    return ControlSet(data)
end

function build_daily_interpolation(time_dict)
    times = Float64[]
    vals = Float64[]
    for (t_str, v) in time_dict
        push!(times, parse(Float64, t_str) * 3600.0)
        push!(vals, Float64(v))
    end
    p = sortperm(times)
    times, vals = times[p], vals[p]
    # Pad to 0-24h
    if times[1] > 0
        times = vcat(0.0, times); vals = vcat(vals[1], vals)
    end
    if last(times) < 86400
        times = vcat(times, 86400.0); vals = vcat(vals, last(vals))
    end
    return ConstantInterpolation(vals, times, extrapolate=true)
end

function get_val(cs::ControlSet, var::Symbol, date::Date, t_seconds::Float64)
    var_data = cs.data[var]
    idx = findfirst(x -> x.first == date, var_data)
    if idx === nothing
        idx = length(var_data) # Fallback to latest defined date
    end
    return var_data[idx].second(t_seconds)
end

# --- 2. CSV Measured Controls ---
struct MeasuredControls
    data::NamedTuple
end

function load_control_interpolations(filename::String)
    df = DataFrame(CSV.File(filename))
    time_abs = df.time .- df.time[1]
    
    nt = (
        U1  = ConstantInterpolation(df.U1, time_abs),
        U2  = ConstantInterpolation(df.U2, time_abs),
        U3  = ConstantInterpolation(df.U3, time_abs),
        U4  = ConstantInterpolation(df.U4, time_abs),
        U5  = ConstantInterpolation(df.U5, time_abs),
        U6  = ConstantInterpolation(df.U6, time_abs),
        U7  = ConstantInterpolation(df.U7, time_abs),
        U8  = ConstantInterpolation(df.U8, time_abs),
        U9  = ConstantInterpolation(df.U9, time_abs),
        U10 = ConstantInterpolation(df.U10, time_abs),
        U11 = ConstantInterpolation(df.U11, time_abs),
        U12 = ConstantInterpolation(df.U12, time_abs),
        Tpipe = ConstantInterpolation(df.Tpipe, time_abs)
    )
    return MeasuredControls(nt)
end

# Helper to access measured controls easily
Base.getproperty(mc::MeasuredControls, s::Symbol) = getproperty(getfield(mc, :data), s)

end # module Controls