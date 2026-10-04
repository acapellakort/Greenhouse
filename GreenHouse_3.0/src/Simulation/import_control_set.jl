using JSON
using DataInterpolations

# Include external utility functions for control
include("../../src/Models/controlFunctions.jl")  

# Set the file path for the controls JSON file
controls_file_path = "src/Simulation/Controls/climateComputerOrders.json"  # Path to your JSON file

# Initialize an empty dictionary to store control instructions
const_dict = Dict()  

# Attempt to read and parse the JSON file into a dictionary
try
    global instruction_dict = JSON.parsefile(controls_file_path)  # Parse the JSON file
catch e
    # Handle any errors that occur during file reading
    println("Error reading JSON file: ", e)
end

# Process the parsed instruction dictionary
for (key, value) in instruction_dict
    try
        println(key)  # Print the current control key
        for (dateString, value_dia) in value
            # Convert the control key to a Symbol and dynamically define the function
            function_name = Symbol(key)
            println(dateString) 
            eval(:($function_name = day_instruction_function($value_dia)))  # Call the day_instruction_function for each day's instructions
        end
    catch e
        # Handle any errors that occur during the processing of control set points
        println("Error processing control set points $key: ", e)
    end
end
########################
