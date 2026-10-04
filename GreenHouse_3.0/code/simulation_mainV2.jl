# In this version rhsPLL_PID! is updated to save all control values
# Probably this code will generate all data and save it as experimental 
# measurements. Th inference is made elsewhere. 


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
# The RHS of greenhouse model
include("utils/greenhouseModel.jl")

using .weatherDataLoader
using .controls
using .psychrometrics
using .GreenhouseModel

# Main simulation routine
#---------------------------------------------------------------------------------------
# Set up of simulations parameters
startDate = Dates.DateTime(2017, 9, 15, 0, 0)
noDays = 50
# Import wetaher file 
pathWeatherFile = "data/weatherData/weatherdata_NL_1988-07-11_2018-07.csv"
# Import control setting points file
pathControlOrdersFile = "controlOrders/climateComputerOrders.json" 
#Obs: Growing period has to be within the available weather data
#Obs: all controls and set points have to be defined in the file at least onces


#---------------------------------------------------------------------------------------
# Preprosesing setup
# Auxiliary time variables and vectors
endDate = startDate + Dates.Day(noDays)
startUnixTime = Dates.datetime2unix(startDate)
timeSpan = 0:300:(noDays*24 * 3600)
days = range(0, stop = noDays, length = noDays+1)
# Time span for each day in seconds (24 hours)
t_span = (0, 24 * 3600.0)
# sampling simulations period in 10 
tsample = 10 # min
# Sampling vector
global teval = (t_span[1]:60*tsample:t_span[2])
# Global variable to store simulation's UnixTime
global unixGlobalTime = startUnixTime

#---------------------------------------------------------------------------------------
# Load weather data between the start and end dates.
Weather = weatherDataLoader.loadInterpolatedWeather(pathWeatherFile, startDate, endDate)

#---------------------------------------------------------------------------------------
# Initila condition('u0') and arrays to store the simulation('u')
global u0 = [18+273.15, 23+273.15, 768.7, 1500.0, 0.0, 0.0, 0.0,]# 0.0, 0.0, 0.0]
global u = u0
"""
  Initial conditions for the simulation

  format: u0 = [ T1, T2, C1, V1, tempIntegralError, co2IntegralError, A, I3, U6, U8]
  units:  [u"K", u"K", u"mg * m^-3", u"Pa", u"K s", u"mg * m^-3 * s", u"1", u"1", u"1", u"1" ]
"""

#---------------------------------------------------------------------------------------
# Simulation start
# Loop over the simulation days: there is an error, so me simulate one day less. 
println("Inicio loop de simulación diario")
@time begin
for day_i in days[1:  noDays]

  #println("Running $day_i ")
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
  
  # Define the ODE problem with current parameters `p`
  #prob = ODEProblem(sistemRHS.simulationODErhs!, u0, t_span, p)
  prob = ODEProblem((du, u, p, t) -> GreenhouseModel.simulationODErhs!(du, u, p, t, Weather, ctrls), u0, t_span, p)

  # Solve the ODE with time evaluation points
  sol = solve(prob, BS3(), saveat = teval)

  # Append the solution and update initial conditions for next day
  global u = hcat(u, sol[1:7,:])
  global unixGlobalTime = vcat(unixGlobalTime, startUnixTime.+ teval .+ day_i*86400)
  global u0 = sol[1:7,end]
end
end

# Convert to hours for plotting
tevalDays = 0.0:600.0:(noDays * 24 * 3600.0)  # every 10 min in seconds

# Convert to DateTime
time_vec = startDate .+ Dates.Second.(tevalDays)
lengthTime = length(time_vec)

# Convert to hours since startDate

plt1=plot(u[1,1:lengthTime].-273.15, lw = 3,label=L"T_{Canopy}")
plt1=plot!(u[2,1:lengthTime].-273.15, lw = 3, label=L"T_{Air}")
plt1=plot!(Weather.Tout(unixGlobalTime).-273.15,lw =2, label =L"T_{out}")
display(plt1)

plt2=plot(time_vec,u[3,1:lengthTime], lw = 3, ylabel=L"\frac{mg}{m^3}", label=L"CO_2")
display(plt2)

RH_sim = psychrometrics.relative_humidity.(u[2,1:lengthTime], u[4,1:lengthTime])
plt3=plot(RH_sim , lw = 3,label=L"RH")
plt3=plot!(Weather.RH(unixGlobalTime),lw =2, label =L"RH_{out}")
display(plt3)

plt4=plot(u[4,1:lengthTime], ylabel=L"Pa", lw = 3,label=L"VP")
plt4=plot!(Weather.VPout(unixGlobalTime),lw =2, label =L"VP_{out}")
display(plt4)

plt5 =  plot(time_vec,u[7,1:lengthTime], lw = 3, ylabel=L"\frac{mg}{ m^3 s}", label=L"Assimilates")
display(plt5)

