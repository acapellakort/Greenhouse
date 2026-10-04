using DataInterpolations

# Auxiliary functions for climate control system

# PI control parameters
# Define the parameters for PI (Proportional-Integral) control:
# T_set = 35.0        # Desired setpoint temperature in °C
# Kp = 0.075          # Proportional gain (response to current error)
# Ki = 0.5            # Integral gain (response to accumulated error)
# decay_factor = 0.1  # Decay factor for leaky integral (slows down the accumulation of old errors)
# dt = 0.1            # Time step for updating the integral (seconds)

# Define the PI control function with leaky integral
function pi_control(V_set, V, integral_error, Kp, Ki, decay_factor, dt)
    # Compute the temperature error
    error = V_set - V
    
    # Update the integral of the error using a leaky integral (decays over time)
    integral_error = integral_error * (1 - decay_factor * dt) + error * dt
    
    # PI control law: a combination of proportional and integral control
    U = Kp * error + Ki * integral_error
    
    # Clamp the control output (heater power) to a value between 0 and 100%
    U = clamp(U, 0, 1)
    
    # Return the control output (U) and the updated integral error
    return U, integral_error
end

# Control for window ventilation using a proportional band (pBand)
# pBand = 8     # Width of proportional band in °C
# ofset = 1     # Offset to avoid premature window opening (°C)
function vent_control(T_set, T2, pBand, ofset )
    # Compare the current temperature (T2) with the setpoint (T_set) + offset

    if T2 < T_set + ofset
        Uvent = 0  # Windows remain closed if the temperature is below the setpoint + offset
    elseif T_set + ofset <= T2 < T_set + pBand + ofset
        # Gradually open windows as temperature rises within the proportional band
        Uvent = (T2 - T_set) / (pBand + ofset)
    elseif T_set + pBand + ofset <= T2
        Uvent = 1  # Fully open windows if the temperature exceeds the proportional band
    end    
    
    # Split ventilation control into side and roof window openings
    U_side = 2 * min(Uvent, 0.5)  # 0 = fully closed, 1 = fully open (side windows)
    U_roof = 2 * max(0, Uvent - 0.5)  # 0 = fully closed, 1 = fully open (roof windows)
    
    # Return control values for side and roof windows
    return U_side, U_roof
end

# Simple algorithm that enable forced ventilation if all windows are open 
# and the temeprature is still high
function Vent_Forced(U8, T_set, T2)
    if U8 == 1 && T2 > T_set 
        U7= (1-T_set/T2)
    else
        U7 = 0
    end    
    return U7
end


function evaluateGreenhouseControls( 
    t :: Float64,
    I5:: Float64,
    ctrls, 
    state::Tuple{Float64, Float64, Float64, Float64}
     ) :: Tuple{Float64, Float64, Float64, Float64, Float64,  Float64, Float64, Float64, Float64, Float64}

    # Unpack orders and set points  
    ToutMax     = ctrls.ToutMax(t)
    Tset        = ctrls.Tset(t)   
    VentpBand   = ctrls.VentpBand(t)
    ofset       = ctrls.ofset(t)
    Light_on    = ctrls.Light_on(t)
    CO2_set     = ctrls.CO2_set(t)

    # Unpack state variables variables
        T1, T2, C1, V1 = state
    # Unpack   
    # Control variables: U1 - U12 sin U10
    # U1, U2, U3, U4, U5, U7, U9, U10, U12 = U1f(t), U2f(t), U3f(t), U4f(t), U5f(t), U7f(t), U9f(t), U10f(t), U12f(t)
    
    U1  = (I5 < ToutMax ? 1.0 : 0.0)        # Thermal screen control   1: open, 0:close
    U2  = 0.0                             # Control ventilador almoadilla
    U3  = 0.0                             # Control sistema enfriamiento  
    U4  = 0.0                             # Heat blower
    U5  = 0.0                             # Sombreado externo
    U6, U8 = vent_control(Tset, T1, VentpBand, ofset) # Ventilation control using side/roof settings
    U7  = 0.0 #Vent_Forced(U8, U11, T1)      # Forced ventilation control
    U9  = 0.0                             # Fog System
    # Light control (on/off)
    U12 = Light_on 
    
    return U1, U2, U3, U4, U5, U6, U7, U8, U9, U12

end

function computePIDTemperature( 
    t :: Float64,
    ctrls, 
    state::Tuple{Float64, Float64, Float64, Float64},
    tempIntegralError:: Float64, 
     ) :: Tuple{Float64, Float64, Float64 }

    # Unpack Temperature Set points  
    Tset = ctrls.Tset(t)   

    # Unpack state variables variables
    T1, T2, C1, V1 = state

    # PI controller for the heating pipe
    U_HeatPipe, tempIntegralError = pi_control(Tset, T1, tempIntegralError, 0.075, 0.5, 0.5, 0.1)
    
    # Heat pipe temperature
    I3 = U_HeatPipe * HEAT_PIPE + (1 - U_HeatPipe) * T2  
    
    # Update PI integral error rate
    tempIntegralErrorRate = Tset - T1        

    return  I3,  tempIntegralError, tempIntegralErrorRate

end

function computePIDCO2( 
    t :: Float64,
    ctrls, 
    state::Tuple{Float64, Float64, Float64, Float64},
    co2IntegralError :: Float64, 
     ) :: Tuple{Float64, Float64, Float64}

    # Unpack orders and set points  
    CO2_set     = ctrls.CO2_set(t)

    # Unpack state variables variables
    T1, T2, C1, V1 = state

    # PI controler CO2
    U10, co2IntegralError = pi_control(CO2_set, C1, co2IntegralError, 0.0030055, 0.00100, 0.5, 0.01)
    
    # Update PI integral error rate
    co2IntegralErrorRate = CO2_set - C1       

    return  U10,  co2IntegralError, co2IntegralErrorRate

end

# This function takes a dictionary of time-value pairs representing
# instructions for a 24-hour period (in seconds) and constructs a 
# step function interpolation over that time.

# Example of input format:
# Dict{String, Any}("6.0" => 0.025966532819438726, "22.0" => 0.13494670391082764, 
#                   "2.0" => 0.02545138681307435, "14.0" => 0.1301136314868927, 
#                   "18.0" => 1.9019747152924538, "10.0" => 0.04047139920294284)

function day_instruction_function(day_instructions)

    # Initialize empty vectors to store time (in seconds) and corresponding values
    values = Float64[]
    times = Float64[]
    
    # Process each time-value pair in the dictionary
    for (time_str, value) in day_instructions
        # Convert time string (hours) to Float64 and multiply to get seconds
        time_hr = parse(Float64, time_str)
        time_sec = time_hr * 60 * 60  # Convert hours to seconds
        
        # Append time in seconds and corresponding value to the vectors
        push!(times, time_sec)
        push!(values, value)
    end
    
    # Sort time and value vectors based on time
    perm_indices = sortperm(times)
    times = times[perm_indices]
    values = values[perm_indices]

    # Ensure the first value starts at 0 seconds (beginning of the day)
    if times[1] != 0.0
        times = vcat(0.0, times)
        values = vcat(values[1], values)
    end
    
    # Ensure the last value extends to the end of the 24-hour period (in seconds)
    if last(times) != 24.0 * 60 * 60
        times = vcat(times, 24.0 * 60 * 60)
        values = vcat(values, last(values))
    end

    # Construct a step function interpolation with constant value between time points
    func = DataInterpolations.ConstantInterpolation(values, times, extrapolate=true)
    
    # Alternatively, we can use linear interpolation if needed:
    # func = DataInterpolations.LinearInterpolation(values, times, extrapolate=true)
    
    # Return the interpolation function
    return func
end


function dayOrders(dateString, instructions, lastOrder)

    # Initialize empty vectors to store time (in seconds) and corresponding values
    orders = Float64[]
    times = Float64[]
    day, month, year = parse.(Int,split(dateString))

    # Process each time-value pair in the dictionary
    for (hourStamp, value) in instructions
        # Convert time string (hours) to Float64 and multiply to get seconds
        hour = parse(Float64, hourStamp)
        dt = Dates.DateTime(year, month, day, hour)
        unixTime = Dates.datetime2unix(dt)
        # Append time in seconds and corresponding value to the vectors
        push!(times, unixTime)
        push!(orders, value)
    end
    
    # Sort time and value vectors based on time
    permIndices = sortperm(times)
    times = times[permIndices]
    orders = orders[permIndices]

    # First hour 
    dt = Dates.DateTime(year, month, day, 0)
    initTime = Dates.datetime2unix(dt)
    # Last hour
    dt = Dates.DateTime(year, month, day, 24)
    lastTime = Dates.datetime2unix(dt)
    # Ensure the first value starts at 0 seconds (beginning of the day)
    if times[1] != initTime
        times = vcat(iniTime, times)
        orders = vcat(lastOrder, orders)
    end
    
    # Ensure the last value extends to the end of the 24-hour period (in seconds)
    if  last(times) != lastTime
        times = vcat(times, lastTime)
        orders = vcat(orders, last(values))
    end

    # Construct a step function interpolation with constant value between time points
    func = DataInterpolations.ConstantInterpolation(values, times, extrapolate = true)
    
    # Alternatively, we can use linear interpolation if needed:
    # func = DataInterpolations.LinearInterpolation(values, times, extrapolate=true)
    
    # Return the interpolation function
    return func
end    
