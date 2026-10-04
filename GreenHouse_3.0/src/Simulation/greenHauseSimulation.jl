# In this version rhsPLL_PID! is updated to save all control values
# Probably this code will generate all data and save it as experimental 
# measurements. Th inference is made elsewhere. 


# Required packages
using CSV
using DataFrames
using Dates
# using Unitful
using DifferentialEquations
using LaTeXStrings

# weather data loader
include("../dataProcessing/weatherData/weatherDataLoader.jl")
# control functions Loader
include("controlLoader.jl")
# Model functions (parameters are defined within these files)
# climate 
include("../Models/climateEquationsRHS.jl")
# photosynthesis
include("../Models/photosynthesisEquations.jl")
# functions for the controls
include("../Models/controlFunctions.jl")
#psychrometrics
# include("../../src/Models/psychrometrics.jl")


using .weatherDataLoader
using .controls

# Define the right-hand side function of the ODE system governing the greenhouse control
# we inlcude the PID controls
function simulationODErhsPID!(du, u, p, t)
  # Unpack variables: [T1, T2, C1, V1, integral_error]
  T1 = u[1]               # Canopy temperature
  T2 = u[2]               # Air temperature
  C1 = u[3]               # CO2 concentration
  V1 = u[4]               # Vapor pressure
  tempIntegralError = u[5]   # Integral error for the PI controller
  co2IntegralError  = u[6]  # Integral error for the PI controller C

  # Convert simulation time `t` to calendar time
  time_absolute = t + p

  # Retrieve weather and environmental data based on the time
  LAI = 0.5 # LAI #I1
  I2 = Weather.Iglob(time_absolute)
  # I3 defined in terms of PI control
  I4 = Weather.Tsky(time_absolute)
  I5 = Weather.Tout(time_absolute)
  I6 = Weather.TmechCool(time_absolute)
  I7 = Weather.Tsoil(time_absolute)
  I8 = Weather.WindSpeed(time_absolute)
  I9 = (1 - eta1) * tau1 * eta2 * Weather.Idocel(time_absolute)
  I10 = Weather.CO2out(time_absolute)
  I11 = Weather.VPout(time_absolute)

  # Photosynthesis assimilation calculation
  A = compute_assimilates(T2, C1, I9, LAI) # Photosynthesis response

 #=
  "Tset"
  "VentpBand"
  "ofset"
  "ToutMax"
  "Light_on"
  "CO2_set"
  "TouMax" # Max temp to open the roof windows
 =#
  
  ToutMax     = ctrls.ToutMax(t)
  Tset        = ctrls.Tset(t)   
  VentpBand   = ctrls.VentpBand(t)
  ofset       = ctrls.ofset(t)
  Light_on    = ctrls.Light_on(t)
  CO2_set     = ctrls.CO2_set(t)
  # Control variables: U1 - U12
  # U1, U2, U3, U4, U5, U7, U9, U10, U12 = U1f(t), U2f(t), U3f(t), U4f(t), U5f(t), U7f(t), U9f(t), U10f(t), U12f(t)
  
    
  U1  = (I5 < ToutMax ? 1.0 : 0.0)        # Thermal screen control   1: open, 0:close
  U2  = 0.0                             # Control ventilador almoadilla
  U3  = 0.0                             # Control sistema enfriamiento  
  U4  = 0.0                             # Heat blower
  U5  = 0.0                             # Sombreado externo
  U6, U8 = vent_control(Tset, T1, VentpBand, ofset) # Ventilation control using side/roof settings
  U7  = 0.0 #Vent_Forced(U8, U11, T1)      # Forced ventilation control
  U9  = 0.0                             # Fog System
  # PI controler CO2
  U10, co2IntegralError = pi_control(CO2_set, C1, co2IntegralError, 0.0030055, 0.00100, 0.5, 0.01)
  # PI controller for the heating pipe
  U_HeatPipe, tempIntegralError = pi_control(Tset, T1, tempIntegralError, 0.075, 0.5, 0.5, 0.1)
  # Heat pipe temperature
  I3 = U_HeatPipe * HEAT_PIPE + (1 - U_HeatPipe) * T2  
  # Light control (on/off)
  U12 = Light_on 

  # Packt State 
  state = (T1, T2, C1, V1)
  
  # Packt controls 
  controls = [U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12]
  
  # Pack inputs
  inputs  = [I2, I3, I4, I5, I6, I7, I8, I9, I10, I11]
      
  # Compute right-hand sides (RHS) of the ODEs
  du[1], du[2], du[3], du[4] = climate_model_rhs(state, A, LAI, inputs, controls)

  # Update integral error and additional variables
  du[5] = Tset - T1        # Integral of error for temperature PI control
  u[5] = tempIntegralError
  du[6] = CO2_set - C1      # Integral of error for CO2 PI control
  u[6] = co2IntegralError   

  # Other variables to keep track of  
  #u[7]  = A                # Photosynthesis response
  #u[8]  = I3               # Heat pipe temperature
  #u[9]  = U6               # Side Ventilation control
  #u[10] = U8               # Roof Ventilation control

end

# Main simulation routine

  # Set up of simulations parameters

  # Simulation´s startdate and days to simulate
  start_date = Dates.DateTime(2018, 8, 15, 10, 0)
  start_time = Dates.datetime2unix(start_date)
  No_days = 5
  days = range(0, stop = No_days, length = No_days + 1)

  # Import weather date. Obs: Growing period has to be within the available weather data
  file_path = "data/weatherdata_NL_1988-07_2018-07.csv"
  end_date = Dates.DateTime(2018, 8, 15 + No_days, 10, 0)
  Weather = weatherDataLoader.loadInterpolatedWeather(file_path, start_date, end_date)

  # Import control setting point functions: Obs: The dates have to be consistent with growing period 
  pathControlOrdersFile = "src/Simulation/Controls/climateComputerOrders.json" 
  currentDay = Dates.datetime2unix(start_date)
  #=
  controls.importControls(currentDay, pathControlOrdersFile)
  ctrls =(
    ToutMax = controls.ToutMax, 
    Tset = controls.Tset,controls.VentpBand,
    ofset = controls.ofset,
    Light_on = controls.Light_on,
    CO2_set = controls.CO2_set)
    =#
  """
    Initial conditions for the simulation

    format: u0 = [ T1, T2, C1, V1, tempIntegralError, co2IntegralError, A, I3, U6, U8]
    units:  [u"K", u"K", u"mg * m^-3", u"Pa", u"K s", u"mg * m^-3 * s", u"1", u"1", u"1", u"1" ]
  """
  global u0 = [18+273.15, 23+273.15, 825.0, 200.0, 0.0, 0.0] #, 0.0, 0.0, 0.0, 0.0]
  
  # Arrays to store the solution and time variable
  global u = u0
  global time_global = start_time
  
  # Time span for each day in seconds (24 hours)
  t_span = (0, 24 * 3600.0)
  # Define time evaluation points every 10 minutes (600 seconds)
  global teval = (t_span[1]:600:t_span[2])

  # Loop over the simulation days
  println("Inicio loop diario")
  for day_i in days
      p = start_time + 86400 * day_i  # Convert day index to Unix time

      println("Running $day_i ")
      # Update day's controls 
      controls.importControls(p, pathControlOrdersFile)
      global ctrls =(
        ToutMax = controls.ToutMax, 
        Tset = controls.Tset,
        VentpBand = controls.VentpBand,
        ofset = controls.ofset,
        Light_on = controls.Light_on,
        CO2_set = controls.CO2_set)
      du = [0.0,0.0,0.0,0.0,0.0,0.0]
      r = simulationODErhsPID!(du, u0, p, t_span[1])
      println(du)
      # Define the ODE problem with current parameters `p`
      prob = ODEProblem(simulationODErhsPID!, u0, t_span, p)
      # Solve the ODE with time evaluation points
      sol = solve(prob, BS3(), saveat = teval)
      #sol = u0
      # Append the solution and update initial conditions for next day
      global u = hcat(u, sol[1:6,:])
      global time_global = vcat(time_global, start_time.+teval .+ day_i*86400)
      global u0 = sol[1:6,end]
  end


