module GreenhouseModel
# Required packages
include("../utils/controlLoader.jl")
# Model functions (parameters are defined within these files)
# climate 
include("../models/climateEquationsRHS.jl")
# photosynthesis
include("../models/photosynthesisEquations.jl")
# functions for the controls
include("../models/controlFunctions.jl")
export simulationODErhs!

  function simulationODErhs!(du, u, p, t, weather, control)
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

    I2 = weather.Iglob(timeAbsolute)
    # I3 defined in terms of PI control
    I4 = weather.Tsky(timeAbsolute)
    I5 = weather.Tout(timeAbsolute)
    I6 = weather.TmechCool(timeAbsolute)
    I7 = weather.Tsoil(timeAbsolute)
    I8 = weather.WindSpeed(timeAbsolute)
    # I9 is the PAR is the par radiation thta reaches the Docel   
    I9 = (1 - eta1) * tau1 * eta2 * weather.Idocel(timeAbsolute) # no lights
    I10 = weather.CO2out(timeAbsolute)
    I11 = weather.VPout(timeAbsolute)


    #---------------------------------------------------------------------------------------
    # Packt State 
    state = (T1, T2, C1, V1)
    # Set U controls, I3 and and PIDerrors
    U1, U2, U3, U4, U5, U6, U7, U8, U9, U12 = evaluateGreenhouseControls(t, I5,control, state)
    #---------------------------------------------------------------------------------------
    # Photosynthesis assimilation calculation
    Itotal = I9 + alpha12*U12
    A = compute_assimilates(T2, C1, Itotal, LAI) # Photosynthesis response
    #--------------------------------------------------------------------------------
    # Set  I3 (pipe temperature) and PID for temperature
    I3,  temperatureIntegralError, temperatureIntegralErrorRate = computePIDTemperature(t, control, state, temperatureIntegralError)
    # Set PID for CO2
    U10, co2IntegralError, co2IntegralErrorRate = computePIDCO2(t, control, state, co2IntegralError)
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


end