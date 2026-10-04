using Plots
using Dates

include("weatherDataLoader.jl")
using .weatherDataLoader

export PlotWeather

"""
xAxisScale(timeSpan::AbstractArray, unit::String="s") -> ScaledTimeAxis

Converts the given time span into an appropriate time scale based on the specified unit. 

### Arguments
- `timeSpan::AbstractArray`: An array of time values to convert.
- `unit::String="s"`: The time unit for scaling. Supported units are:
  - "s": Seconds
  - "min": Minutes
  - "h": Hours

### Returns
- `ScaledTimeAxis`: An array of time values converted to the specified time unit.

### Example
julia> xAxisScale(0:300:7200, "h")
# Returns an array of values representing time in hours, based on the given `timeSpan`
"""
function xAxisScale(timeSpan::AbstractArray, unit::String="s")
    conversionFactor = Dict("s" => 1, "min" => 60, "h" => 3600)
    return timeSpan ./ conversionFactor[unit]
end


"""
PlotWeather(keys::Vector{Symbol}, startDate::DateTime, timeSpan::AbstractArray, weather, 
            timeUnit::String="min")

Generates plots for specified weather variables over a given time span.

### Arguments
- `keys::Vector{Symbol}`: A vector of symbols representing the weather data keys to plot.
- `startDate::DateTime`: The starting date and time for the data.
- `timeSpan::AbstractArray`: An array representing the time span (in seconds) for which the data will be plotted.
- `weather::WeatherData`: An instance of `WeatherData` that contains the necessary weather data.
- `timeUnit::String="min"`: The unit to use for scaling the time axis. Supported values are:
  - "min": Minutes
  - "h": Hours
  - "s": Seconds (default)

### Keyword Arguments
- `timeUnit::String="min"`: Specifies the unit of time for scaling the x-axis. Defaults to "min" (minutes).

### Example
julia> PlotWeather([:Tout, :RH, :Iglob, :Tpipe], start_date, timeSpan, Weather, timeUnit="h")

This function generates a plot for each weather variable specified in `keys`, showing the data over the given `timeSpan` in the specified `timeUnit`. The plot displays time on the x-axis and the corresponding weather variable values on the y-axis.

### Notes
- The `key` parameter determines the specific type of weather data to plot (e.g., :Tout for temperature, :RH for relative humidity).
- Time is converted to the specified unit using the `xAxisScale` function.
- If the first character of `key` is 'T', the plot will adjust the y-axis to account for temperature in Kelvin by subtracting 273.15.
"""

function PlotWeather(keys::Vector{Symbol}, startDate::DateTime, timeSpan::AbstractArray, weather; timeUnit::String="min")
    t0 = Dates.datetime2unix(startDate)
    for key in keys
        variableFunction = getproperty(weather, key)
        keyName = string(:Tout)
        values = variableFunction.(t0 .+ timeSpan)
        timeAxis = xAxisScale(timeSpan, timeUnit)
        
        if string(key)[1] == 'T'
            yShift = 273.15
        else
            yShift = 0
        end
        
        p = plot(timeAxis, values .- yShift, xlabel="time", label=string(key))
        display(p)
    end
end

file_path = joinpath(@__DIR__, "..", "data", "weatherData", "weatherdata_NL_1988-07-11_2018-07.csv")

#file_path = "../data/weatherData/weatherdata_NL_1988-07_2018-07.csv"
#=
start_date = Dates.DateTime(2018, 8, 15, 10, 0)
No_days = 10
end_date = Dates.DateTime(2018, 8, 15 + No_days, 10, 0)
Weather = weatherDataLoader.loadInterpolatedWeather(file_path, start_date, end_date)
timeSpan = 0:300:(No_days*24 * 3600)  # every 600 segs (10 min) for one day

keys = [ :Tout, :RH, :Iglob, :Tpipe ]

#PlotWeather(keys, start_date, timeSpan, Weather, timeUnit="h")
=#