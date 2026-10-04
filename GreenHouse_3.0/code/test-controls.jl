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
startDate = Dates.DateTime(2018, 8, 15, 0, 0)
noDays = 5
endDate = startDate + Dates.Day(noDays)
startUnixTime = Dates.datetime2unix(startDate)
timeSpan = 0:300:(noDays*24 * 3600)
days = range(0, stop = noDays, length = noDays+1)
#---------------------------------------------------------------------------------------
# Import control setting point functions: Obs: The dates have to be consistent with growing period 
pathControlOrdersFile = "controlOrders/climateComputerOrders.json" 
  
for day_i in days[1:  noDays]

    println("Running $day_i ")
    p = startUnixTime + 86400 * day_i  # Update starting daily simulation date 

    # Update day's controls 
    controls.importControls(p, pathControlOrdersFile)
    ctrls =(
    ToutMax = controls.ToutMax, 
    Tset = controls.Tset,
    VentpBand = controls.VentpBand,
    ofset = controls.ofset,
    Light_on = controls.Light_on,
    CO2_set = controls.CO2_set)

    plt = plot(ctrls.Tset)
    display(plt) 
end
