using Dates
using JSON

function json_to_vectors(json_data::String)
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
            datetime = DateTime(year(date), month(date), day(date), hour, minute)
            
            # Convert DateTime to Unix time
            unix_time = Dates.datetime2unix(datetime)
            
            # Append to vectors
            push!(unix_times, unix_time)
            push!(values, value)
        end
    end

    # Get permutation indices to sort unix_times
    perm_indices = sortperm(unix_times)

    # Apply permutation to both vectors
    unix_times = unix_times[perm_indices]
    values = values[perm_indices]

    return unix_times, values
end

using DataInterpolations

function create_interpolation_function(unix_times::Vector{Float64}, values::Vector{Float64})
    # Create an interpolation object
    interpolant = DataInterpolations.ConstantInterpolation( values, unix_times)

    # Return the interpolation function
    return interpolant
end

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
            "10.0": 1.028020287514664236,
            "14.0": 0.020747496630065147,
            "18.0": 1.5197467480320483,
            "22.0": 0.004685055144364013
        }
    }
}
"""

# Convert JSON data to vectors
unix_times, values = json_to_vectors(json_data)

# Create interpolation function
interpolant = create_interpolation_function(unix_times, values)

# Use the interpolation function
sample_time = Dates.datetime2unix(DateTime(2024, 9, 10, 12, 30))
interpolated_value = interpolant(sample_time)

println("Interpolated value at $sample_time: $interpolated_value")


