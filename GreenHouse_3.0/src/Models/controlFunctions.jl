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
