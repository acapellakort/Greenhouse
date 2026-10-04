
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
    itp = DataInterpolations.ConstantInterpolation(orders, times; 
            extrapolation_left = DataInterpolations.ExtrapolationType.Constant,
            extrapolation_right = DataInterpolations.ExtrapolationType.Constant
        )
    #func = DataInterpolations.ConstantInterpolation(orders, times, extrapolate = true)

    # Alternatively, we can use linear interpolation if needed:
    # func = DataInterpolations.LinearInterpolation(values, times, extrapolate=true)
    
    # Return the interpolation function
    return itp
end    


"""
    find_date_range(dict::Dict{String, Any}) -> (Date, Date)

Given a nested dictionary structure where the second-level keys are date strings
in the "dd-mm-yyyy" format, this function extracts all such date strings, parses 
them into `Date` objects, and returns a tuple containing the earliest and latest dates.

# Arguments
- `dict::Dict{String, Any}`: A dictionary where some values are nested dictionaries 
  keyed by date strings.

# Returns
- A tuple `(earliest_date::Date, latest_date::Date)` representing the range of dates found.
"""
function find_date_range(dict::Dict{String, Any})
    dates = String[]
    
    for (key, subdict) in dict
        if isa(subdict, Dict)
            append!(dates, keys(subdict))
        end
    end

    # Remove duplicates and parse dates
    unique_dates = unique(dates)
    parsed_dates = Date.(unique_dates, dateformat"dd-mm-yyyy")

    earliest = Dates.format(minimum(parsed_dates), dateformat"dd-mm-yyyy")
    latest = Dates.format(maximum(parsed_dates), dateformat"dd-mm-yyyy")
    return earliest, latest
end

"""
    find_date_range_per_key(dict::Dict{String, Any}) -> Dict{String, Tuple{String, String}}

For each top-level key in the input dictionary, finds the earliest and latest date strings
(keyed as "dd-mm-yyyy") from its nested dictionary.

# Arguments
- `dict::Dict{String, Any}`: A dictionary where each value is a nested dictionary with date strings as keys.

# Returns
- A dictionary mapping each top-level key to a tuple of strings `(earliest_date, latest_date)`,
  both in "dd-mm-yyyy" format.
"""
function find_date_range_per_key(dict::Dict{String, Any})
    result = Dict{String, Tuple{String, String}}()

    for (outer_key, inner_dict) in dict
        if isa(inner_dict, Dict)
            date_keys = collect(keys(inner_dict))
            parsed_dates = Date.(date_keys, dateformat"dd-mm-yyyy")
            earliest = Dates.format(minimum(parsed_dates), dateformat"dd-mm-yyyy")
            latest = Dates.format(maximum(parsed_dates), dateformat"dd-mm-yyyy")
            result[outer_key] = (earliest, latest)
        end
    end

    return result
end


"""
    get_date_list_per_key(dict::Dict{String, Any}) -> Dict{String, Vector{String}}

Returns a dictionary mapping each top-level key (e.g., "Tset") to a sorted list of date strings 
in "dd-mm-yyyy" format that are found in its nested dictionary.
"""
function get_date_list_per_key(dict::Dict{String, Any})
    result = Dict{String, Vector{String}}()

    for (outer_key, inner_dict) in dict
        if isa(inner_dict, Dict)
            dates = sort(collect(keys(inner_dict)))  # keep them sorted as strings
            result[outer_key] = dates
        end
    end

    return result
end

"""
    select_date(current_date_str::String, date_list::Vector{String}) -> String

Selects the most recent date from `date_list` (in "dd-mm-yyyy" format) that is less than or equal to
`current_date_str`. If none are earlier or equal, returns the first date in the list.
"""
function select_date(current_date_str::String, date_list::Vector{String})
    current_date = Date(current_date_str, dateformat"dd-mm-yyyy")
    parsed_dates = sort(Date.(date_list, dateformat"dd-mm-yyyy"))

    # Filter dates less than or equal to current_date
    candidates = filter(d -> d <= current_date, parsed_dates)

    selected = isempty(candidates) ? parsed_dates[1] : candidates[end]
    return Dates.format(selected, dateformat"dd-mm-yyyy")
end

#---------------------------------------------------------------------------------

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
    # Recover date format for current day
    currentDayKey = Dates.format(Dates.unix2datetime(currentDay),"dd-mm-yyyy")
    # Find the earliest and latest instriction day
    controlDatesPerKey = get_date_list_per_key(ordersDict) 


    # Process the parsed instruction dictionary
    for (controlName, instructions) in ordersDict
        try 
            controlSelectedDate = select_date(currentDayKey, controlDatesPerKey[controlName] )
            dayInstructions = instructions[controlSelectedDate]
            
            # Call the dayOrders for controlSelectedDate's instructions
            eval(:(global $(Symbol(controlName)) = dayOrders($dayInstructions)))

        catch e
            # Handle any errors that occur during the processing of control set points
            println("Error processing control set points $controlName: ", e)
        end   
    end
end


function usage()
    pathOrdersFile = "controlOrders/climateComputerOrders.json" 
    currentDay = Dates.datetime2unix(Dates.DateTime(2018,7,15))
    ctrl = importControls(currentDay, pathOrdersFile)

end
end 

#using .controls
#pathOrdersFile = "controlOrders/climateComputerOrders.json" 
#currentDay = Dates.datetime2unix(Dates.DateTime(2024,09,11))
#controls.importControls(currentDay, pathOrdersFile)

#plot(controls.Tset)

#=

d = Dict{String, Any}(...)  # Your dictionary here
earliest, latest = find_date_range(d)
println("Earliest: $earliest, Latest: $latest")

=#