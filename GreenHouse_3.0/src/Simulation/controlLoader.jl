
module controls
using Dates
using JSON
using DataInterpolations


export importControls


"""
dayInstruction(instructions) 

Instructions is a dictionary of hourString => value 

"""
function dayOrders(instructions)

    # Initialize empty vectors to store time (in seconds) and corresponding values
    orders = Float64[]
    times = Float64[]
   
    # Process each time-value pair in the dictionary
    for (hourStamp, value) in instructions
        # Convert time string (hours) to Float64 and multiply to get seconds

        hour = parse(Float64, hourStamp)
        seconds = 3600*hour
        # Append time in seconds and corresponding value to the vectors
        push!(times, seconds)
        push!(orders, value)
    end
   
    # Sort time and value vectors based on time
    permIndices = sortperm(times)
    times = times[permIndices]
    orders = orders[permIndices]
    

    # Ensure the first value starts at 0 seconds (beginning of the day)
    if times[1] != 0.0
        times = vcat(0.0, times)
        orders = vcat(orders[1], orders)
    end

    # Ensure the last value extends to the end of the 24-hour period (in seconds)
    if  last(times) != 24*3600
        times = vcat(times, 24*3600)
        orders = vcat(orders, last(orders))
    end

    # Construct a step function interpolation with constant value between time points
    func = DataInterpolations.ConstantInterpolation(orders, times, extrapolate = true)

    # Alternatively, we can use linear interpolation if needed:
    # func = DataInterpolations.LinearInterpolation(values, times, extrapolate=true)
    
    # Return the interpolation function
    return func
end    



######################

function loadDict(pathOrdersFile)
    try
        # Parse the JSON file as dictionary to store orders
        ordersDict = JSON.parsefile(pathOrdersFile)   
        return ordersDict
    catch e
        # Handle any errors that occur during file reading
        println("Error reading JSON file: ", e)
    end
    

end 

function importControls(currentDay, pathOrdersFile)

# Attempt to read and parse the JSON file into a dictionary
ordersDict = loadDict(pathOrdersFile)

# Process the parsed instruction dictionary
for (controlName, instructions) in ordersDict
    #functionName = Symbol(controlName)
    try 
        currentDayKey = Dates.format(Dates.unix2datetime(currentDay),"dd-mm-yyyy")
        if haskey(instructions, currentDayKey)
            dayInstructions = instructions[currentDayKey]
            # Call the dayOrders for currentday's instructions
            eval(:(global $(Symbol(controlName)) = dayOrders($dayInstructions)))
        end
       
    catch e
        # Handle any errors that occur during the processing of control set points
        println("Error processing control set points $controlName: ", e)
    end   

end


end


function usage()
    pathOrdersFile = "src/Simulation/Controls/climateComputerOrders.json" 
    currentDay = Dates.datetime2unix(Dates.DateTime(2024,09,11))
    importControls(currentDay, pathOrdersFile)
end
end 

#using .controls
#pathOrdersFile = "src/Simulation/Controls/climateComputerOrders.json" 
#currentDay = Dates.datetime2unix(Dates.DateTime(2024,09,11))
#controls.importControls(currentDay, pathOrdersFile)

#plot(controls.Tset)