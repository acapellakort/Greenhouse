# Required packages
using CSV
using DataFrames
using Dates
using Plots
# using Unitful
using DifferentialEquations
using LaTeXStrings


# weather data loader
include("utils/weatherDataLoader.jl")
# control functions Loader
include("utils/controlLoader.jl")
# Model functions (parameters are defined within these files)
# climate 
include("models/climateEquationsRHS.jl")
# photosynthesis
include("models/photosynthesisEquations.jl")
# functions for the controls
include("models/controlFunctions.jl")
#psychrometrics
# include("../../src/Models/psychrometrics.jl")

include("utils/plotUtils.jl")

using .weatherDataLoader
using .controls


#---------------------------------------------------------------------------------------
# Simulation's startdate and days to simulate
start_date = Dates.DateTime(2017, 8, 15, 0, 0)
No_days = 12
end_date = start_date + Dates.Day(No_days)
timeSpan = 0:300:(No_days*24 * 3600)

#---------------------------------------------------------------------------------------
# Import weather. Obs: Growing period has to be within the available weather data
pathWeatherFile = "data/weatherData/weatherdata_NL_1988-07-11_2018-07.csv"

Weather = weatherDataLoader.loadInterpolatedWeather(pathWeatherFile, start_date, end_date)

keys = [ :Tout, :RH, :Iglob, :VPout]
PlotWeather(keys, start_date, timeSpan, Weather, timeUnit="h")