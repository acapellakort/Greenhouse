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


using .weatherDataLoader
using .controls
using .psychrometrics


# Define the right-hand side function of the ODE system governing the greenhouse control
# we inlcude the PID controls
function simulationODErhsPID!(du, u, p, t)
  # Unpack variables: [T1, T2, C1, V1, integral_error]
  T1 = u[1]                         # Canopy temperature
  T2 = u[2]                         # Air temperature
  C1 = u[3]                         # CO2 concentration
  V1 = u[4]                         # Vapor pressure
  temperatureIntegralError = u[5]   # Integral error for the PI controller
  co2IntegralError  = u[6]          # Integral error for the PI controller C

  # Convert simulation time `t` to calendar time
  timeAbsolute = t + p

  # Retrieve weather and environmental data based on the time
  LAI = 0.5 # LAI #I1

  I2 = Weather.Iglob(timeAbsolute)
  # I3 defined in terms of PI control
  I4 = Weather.Tsky(timeAbsolute)
  I5 = Weather.Tout(timeAbsolute)
  I6 = Weather.TmechCool(timeAbsolute)
  I7 = Weather.Tsoil(timeAbsolute)
  I8 = Weather.WindSpeed(timeAbsolute)
  # I9 is the PAR is the par radiation thta reaches the Docel   
  I9 = (1 - eta1) * tau1 * eta2 * Weather.Idocel(timeAbsolute) # No ligths
  I10 = Weather.CO2out(timeAbsolute)
  I11 = Weather.VPout(timeAbsolute)

  
  #---------------------------------------------------------------------------------------
  # Packt State 
  state = (T1, T2, C1, V1)
  # Set U controls, I3 and and PIDerrors
  U1, U2, U3, U4, U5, U6, U7, U8, U9, U12 = evaluateGreenhouseControls(t, I5,ctrls, state)
  #---------------------------------------------------------------------------------------
  # Photosynthesis assimilation calculation
  Itotal = I9 + alpha12*U12
  A = compute_assimilates(T2, C1, Itotal, LAI) # Photosynthesis response
  #---------------------------------------------------------------------------------------
  # Set  I3 (pipe temperature) and PID for temperature
  I3,  temperatureIntegralError, temperatureIntegralErrorRate = computePIDTemperature(t, ctrls, state, temperatureIntegralError)
  # Set PID for CO2
  U10, co2IntegralError, co2IntegralErrorRate = computePIDCO2(t, ctrls, state, co2IntegralError)
  # Update integral error and additional variables           
  u[5]  = temperatureIntegralError               # Integral of error for temperature PI control
  du[5] = temperatureIntegralErrorRate    # Integral of error Rate for temperature PI control
  u[6]  = co2IntegralError                # Integral of error for CO2 PI control
  du[6] = co2IntegralErrorRate            # Integral of error Rate for CO2 PI control
  #---------------------------------------------------------------------------------------
  # Packt controls 
  controls = [U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12]
  # Pack inputs
  inputs  = [I2, I3, I4, I5, I6, I7, I8, I9, I10, I11]
  #---------------------------------------------------------------------------------------
  # Compute right-hand sides (RHS) of the ODEs
  du[1], du[2], du[3], du[4] = climate_model_rhs(state, A, LAI, inputs, controls)

  # Other variables to keep track of  
  u[7]  = A                # Photosynthesis response
  #u[8]  = I3               # Heat pipe temperature
  #u[9]  = U6               # Side Ventilation control
  #u[10] = U8               # Roof Ventilation control

end

# Main simulation routine

  # Set up of simulations parameters

  #---------------------------------------------------------------------------------------
  # Simulation´s startdate and days to simulate
  startDate = Dates.DateTime(2017, 9, 15, 0, 0)
  noDays = 90
  endDate = startDate + Dates.Day(noDays)
  startUnixTime = Dates.datetime2unix(startDate)
  timeSpan = 0:300:(noDays*24 * 3600)
  days = range(0, stop = noDays, length = noDays+1)
  #---------------------------------------------------------------------------------------
  # Import weather. Obs: Growing period has to be within the available weather data
  pathWeatherFile = "data/weatherData/weatherdata_NL_1988-07-11_2018-07.csv"
  
  Weather = weatherDataLoader.loadInterpolatedWeather(pathWeatherFile, startDate, endDate)

  #---------------------------------------------------------------------------------------
  # Import control setting point functions: Obs: The dates have to be consistent with growing period 
  pathControlOrdersFile = "controlOrders/climateComputerOrders.json" 
  

  #---------------------------------------------------------------------------------------
  """
    Initial conditions for the simulation

    format: u0 = [ T1, T2, C1, V1, tempIntegralError, co2IntegralError, A, I3, U6, U8]
    units:  [u"K", u"K", u"mg * m^-3", u"Pa", u"K s", u"mg * m^-3 * s", u"1", u"1", u"1", u"1" ]
  """
  global u0 = [18+273.15, 23+273.15, 768.7, 1500.0, 0.0, 0.0, 0.0]#, 0.0, 0.0, 0.0]
  
  #---------------------------------------------------------------------------------------
  # Arrays to store the solution and time variable
  global u = u0
  global unixGlobalTime = startUnixTime
  
  #---------------------------------------------------------------------------------------
  # Time span for each day in seconds (24 hours)
  t_span = (0, 24 * 3600.0)
  # Define time evaluation points every 10 minutes (600 seconds)
  global teval = (t_span[1]:600:t_span[2])

  #---------------------------------------------------------------------------------------
  # Simulation start
  # Loop over the simulation days
  println("Inicio loop diario")
  for day_i in days[1:  noDays]

    println("Running $day_i ")
    p = startUnixTime + 86400 * day_i  # Update starting daily simulation date 


    # Update day's controls 
    controls.importControls(p, pathControlOrdersFile)
    global ctrls =(
    ToutMax = controls.ToutMax, 
    Tset = controls.Tset,
    VentpBand = controls.VentpBand,
    ofset = controls.ofset,
    Light_on = controls.Light_on,
    CO2_set = controls.CO2_set)
    
    # Define the ODE problem with current parameters `p`
    prob = ODEProblem(simulationODErhsPID!, u0, t_span, p)

 
    # Solve the ODE with time evaluation points
    sol = solve(prob, BS3(), saveat = teval)

    # Append the solution and update initial conditions for next day
    global u = hcat(u, sol[1:7,:])
    global unixGlobalTime = vcat(unixGlobalTime, startUnixTime.+ teval .+ day_i*86400)
    global u0 = sol[1:7,end]
  end


# Convert to hours for plotting
tevalDays = 0.0:600.0:(noDays * 24 * 3600.0)  # every 10 min in seconds

# Convert to DateTime
time_vec = startDate .+ Dates.Second.(tevalDays)
lengthTime = length(time_vec)

# Convert to hours since startDate

plt=plot(u[1,1:lengthTime].-273.15, lw = 3,label=L"T_{Canopy}")
plt=plot!(u[2,1:lengthTime].-273.15, lw = 3, label=L"T_{Air}")
plt=plot!(Weather.Tout(unixGlobalTime).-273.15,lw =2, label =L"T_{out}")
display(plt)

plt=plot(time_vec,u[3,1:lengthTime], lw = 3, ylabel=L"\frac{mg}{m^3}", label=L"CO_2")
display(plt)

RH_sim = psychrometrics.relative_humidity.(u[2,1:lengthTime], u[4,1:lengthTime])
plt=plot(RH_sim , lw = 3,label=L"RH")
plt=plot!(Weather.RH(unixGlobalTime),lw =2, label =L"RH_{out}")
display(plt)

plt=plot(u[4,1:lengthTime], ylabel=L"Pa", lw = 3,label=L"VP")
plt=plot!(Weather.VPout(unixGlobalTime),lw =2, label =L"VP_{out}")
display(plt)