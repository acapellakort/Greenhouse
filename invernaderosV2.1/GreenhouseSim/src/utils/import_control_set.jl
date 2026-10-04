using JSON
using DataInterpolations

include("../utils/control_functions.jl")  # Include external utility functions for control

base_dir = pwd() 

# Build paths using joinpath (this works on Mac, Linux, and Windows automatically)
controls_file_path = joinpath(base_dir, "configfiles", "control_instructions.json")

# Set the file path for the controls JSON file
#controls_file_path = "/Users/capella/Documents/trabajo/academia/unam/proyectos investigacion/bio/invernaderos/invernaderos2.0/code/configfiles/control_instructions.json"  # Path to your JSON file

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
        for (key2, value_dia) in value
            # Convert the control key to a Symbol and dynamically define the function
            function_name = Symbol(key)
            eval(:($function_name = day_instruction_function($value_dia)))  # Call the day_instruction_function for each day's instructions
        end
    catch e
        # Handle any errors that occur during the processing of control set points
        println("Error processing control set points $key: ", e)
    end
end
