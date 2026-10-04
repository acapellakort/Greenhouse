module weatherDataLoader

using CSV
using DataFrames
using Dates
using Unitful
using DataInterpolations
using Plots

include("../../Models/psychrometrics.jl")
using .psychrometrics

export loadInterpolatedWeather

"""
# weatherDataLoader

This module provides tools to load, filter, and process weather data from CSV files and interpolate the data into functions for simulation purposes.

## Features
- Load weather data from a CSV file.
- Filter weather data based on a specified date range.
- Compute derived columns such as Unix time and vapor pressure.
- Interpolate weather parameters for modeling and analysis.

## Example Usage
```julia
using weatherDataLoader

file_path = "data/weatherdata_NL_1988-07_2018-07.csv"
start_date = Dates.DateTime(2018, 8, 15, 10, 0)
end_date = Dates.DateTime(2018, 8, 17, 10, 0)

Weather = weatherDataLoader.load_interpolated_weather(file_path, start_date, end_date)
plot(Weather.Tout)
```

## Exports
- `load_interpolated_weather`: Load and interpolate weather data for modeling.
- `usage`: Demonstrates example usage of the module.

## Dependencies
- `CSV`: For reading CSV files.
- `DataFrames`: For handling tabular data.
- `Dates`: For working with date and time.
- `Unitful`: For handling data with units.
- `DataInterpolations`: For creating interpolation functions.
- `Plots`: For visualization of data.

"""

# Function to compute Unix time and replace the column
function add_unix_time!(df)
    for row in eachrow(df)
        dt = Dates.DateTime(row.Year, row.Month, row.Day, row.Hour, row.Minute)
        row.Unix_Time = Dates.datetime2unix(dt)
    end
end

# Function to compute vapor pressure and replace the column
function add_Vapor_Pressure!(df)
    for row in eachrow(df)
        VP = psychrometrics.vapor_pressure(row.Temperature, row.Relative_Humidity)
        row.Vapor_Pressure = VP
    end
end

# Function to load and filter weather data
function load_and_filter_weather_data(file_path::String, start_date, end_date)
    df = CSV.read(file_path, DataFrame)

    rename!(df, 1  => "Unix_Time")
    rename!(df, 7  => "Temperature")
    rename!(df, 8  => "Relative_Humidity")
    rename!(df, 9  => "Vapor_Pressure")
    rename!(df, 13 => "External_Radiation")
    rename!(df, 14 => "Wind_Speed")

    add_unix_time!(df)
    add_Vapor_Pressure!(df)

    unix_start_date = Dates.datetime2unix(start_date)
    unix_end_date = Dates.datetime2unix(end_date)
    filtered_df = df[(df.Unix_Time .>= unix_start_date) .& (df.Unix_Time .<= unix_end_date), :]

    keep_columns = ["Unix_Time", "Year", "Month", "Day", "Hour", "Minute", "Temperature", 
                    "Relative_Humidity", "Vapor_Pressure", "External_Radiation", "Wind_Speed"]
    filtered_df = select(filtered_df, keep_columns)

    filtered_df.Unix_Time = filtered_df.Unix_Time
    filtered_df.Temperature = filtered_df.Temperature .+ 273.15
    filtered_df.Vapor_Pressure = filtered_df.Vapor_Pressure
    filtered_df.External_Radiation = filtered_df.External_Radiation
    filtered_df.Wind_Speed = filtered_df.Wind_Speed .* (5.0 / 18.0)

    return filtered_df
end

# Function to load and interpolate weather data
function loadInterpolatedWeather(file_path::String, start_date, end_date)
    filtered_data = load_and_filter_weather_data(file_path, start_date, end_date)

    Toutf = DataInterpolations.ConstantInterpolation(filtered_data.Temperature, filtered_data.Unix_Time)
    RHf = DataInterpolations.ConstantInterpolation(filtered_data.Relative_Humidity, filtered_data.Unix_Time)
    VPoutf = DataInterpolations.ConstantInterpolation(filtered_data.Vapor_Pressure, filtered_data.Unix_Time)
    Iglobalf = DataInterpolations.ConstantInterpolation(filtered_data.External_Radiation, filtered_data.Unix_Time)
    WindSpeedf = DataInterpolations.ConstantInterpolation(filtered_data.Wind_Speed, filtered_data.Unix_Time)

    CO2outf(t) = 834.7
    Tpipef(f) = 15.0 + 273.15
    Tskyf(t) = -0.4 + 273.15
    TmechCoolf(t) = 20.0 + 273.15
    Tsoilf(t) = 18.0 + 273.15

    Idocelf = Iglobalf

    result = (Iglob = Iglobalf, Tpipe = Tpipef, Tsky = Tskyf,
             Tout = Toutf, TmechCool = TmechCoolf, Tsoil = Tsoilf,
             WindSpeed = WindSpeedf, Idocel = Idocelf, CO2out = CO2outf, 
             VPout = VPoutf, RH = RHf)
    return result
end

# Example usage demonstration
function usage()
    println("Example usage:")
    println("file_path = \"data/weatherdata_NL_1988-07_2018-07.csv\"")
    println("start_date = Dates.DateTime(2018, 8, 15, 10, 0)")  
    println("No_days = 2")
    println("end_date = Dates.DateTime(2018, 8, 15 + No_days, 10, 0)")
    println("Weather = weatherDataLoader.loadInterpolatedWeather(file_path, start_date, end_date)")
    println("")
    println("plot(Weather.RH)")

    file_path = "data/weatherdata_NL_1988-07_2018-07.csv"
    start_date = Dates.DateTime(2018, 8, 15, 10, 0)
    No_days = 2
    end_date = Dates.DateTime(2018, 8, 15 + No_days, 10, 0)

    Weather = loadInterpolatedWeather(file_path, start_date, end_date)
    Plots.plot(Weather.Tout)
end

end