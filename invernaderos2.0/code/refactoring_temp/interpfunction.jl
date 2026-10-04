using Dates
using JSON
using Interpolations
using BenchmarkTools

# Sample JSON data
json_data = """
{
    "@minTemp": {
        "10-09-2024": {
            "2.0": 0.02545138681307435,
            "6.0": 0.025966532819438726,
            "10.0": 0.04047139920294284,
            "14.0": 0.1301136314868927,
            "18.0": 1.9019747152924538,
            "22.0": 0.13494670391082764
        },
        "11-09-2024": {
            "2.0": 1.06372811365872622,
            "6.0": 0.010586723918095233,
            "10.0": 0.028020287514664236,
            "14.0": 0.020747496630065147,
            "18.0": 1.5197467480320483,
            "22.0": 0.004685055144364013
        }
    }
}
"""

# Parse the JSON data
parsed_data = JSON.parse(json_data)
min_temp_data = parsed_data["@minTemp"]

# Convert date format
date_format = "dd-mm-yyyy"

# Prepare vectors
unix_times = Float64[]
values = Float64[]

# Process each date in the JSON data
for (date_str, times) in min_temp_data
    # Convert date string to Date object
    date = Date(date_str, date_format)
    
    # Process each time entry
    for (time_str, value) in times
        # Convert time string to Float64 (assuming it is in hours)
        time_float = parse(Float64, time_str)
        
        # Extract hours and minutes
        hour = floor(Int, time_float)
        minute = round(Int, (time_float - hour) * 60)
        
        # Create DateTime object
        datetime2 = DateTime(year(date), month(date), day(date), hour, minute)
        
        # Convert DateTime to Unix time
        unix_time = Dates.datetime2unix(datetime2)
        
        # Append to vectors
        push!(unix_times, unix_time)
        push!(values, value)
    end
end


# Get permutation indices to sort vector1
perm_indices = sortperm(unix_times)

# Apply permutation to both vectors
unix_times = unix_times[perm_indices]
values = values[perm_indices]

# Print the sorted vectors

# Display the vectors
println("Unix Times: ", unix_times)
println("Values: ", values)

using DataInterpolations,  Plots
gr()

unix_times=unix_times

@btime A = DataInterpolations.ConstantInterpolation(values, unix_times, extrapolate=true)
A = DataInterpolations.ConstantInterpolation(values, unix_times, extrapolate=true)

plot(A)


for (date, inner_dict) in data
    # Sort the keys lexicographically
    keys_sorted = sort(collect(keys(inner_dict)))
    println(keys_sorted)
    end



    for (date, inner_dict) in data
        # Sort the keys lexicographically
        keys_sorted = sort(collect(keys(inner_dict)))

        # Check if the first key (lexicographically) is not "0.0"
        if keys_sorted[1] != "0.0"
            # Get the value of the first key
            first_value = inner_dict[keys_sorted[1]]
            # Add a new entry for "0.0" with the first value
            inner_dict["0.0"] = first_value
        end
        
        # Re-sort the dictionary after the insertion
        sorted_dict = Dict(k => inner_dict[k] for k in sort(collect(keys(inner_dict))))
        
        # Update the main dictionary
        data[date] = sorted_dict
    end

    # Print the updated dictionary


    # Process each date in the JSON data
function day_interpolation(day_data)
    # Prepare data
    values = Float64[]
    times = Float64[]
    # Process each time entry
    for (time_str, value) in day_data
        # Convert time string to Float64 (assuming it is in hours)
        time_hr = parse(Float64, time_str)
        
        # Extract time in seconds
        time_sec = time_hr*60*60 
        
        # Append to vectors
        push!(times, time_sec)
        push!(values, value)

    end
    # Get permutation indices to sort vector1
    perm_indices = sortperm(times)

    # Apply permutation to both vectors
    times = times[perm_indices]
    values = values[perm_indices]


    # Complete firsty value if necessary 
    if times[1] != 0.0
        times = vcat(0.0, times)
        values = vcat(values[1], values)
    end
    # complete last value if necessary
    if last(times) != 24.0*60*60
        times = vcat(times, 24.0*60*60)
        values = vcat(values, last(values))
    end

    func = DataInterpolations.ConstantInterpolation(values, times, extrapolate=true)
    return func
end

