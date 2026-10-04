using JSON
using Unitful


# Defining new units
Unitful.register(@__MODULE__);
@unit USD "USD" Dollar 1.0 false
@unit MXN "MXN" Pesos 1.0 false
@unit EUR "EUR" Euros 1.0 false
@unit ppm "ppm" part_per_million 1e-6 false



# Define a mapping from unit strings in your JSON file to Unitful units
unit_mapping = Dict(
    "1" => u"NoUnits",
    "m" => u"m",
    "meters" => u"m",
    "metros" => u"m",
    "m * s^-2" => u"m * s^-2",
    "m * W^-2" => u"m * W^-2",
    "m * m^-2" => u"m * m^-2",
    "m^2" => u"m^2",
    "m^2*m^-2" => u"m^2*m^-2",
    "m^3 * s^-1" => u"m^3 * s^-1",
    "seconds" => u"s",
    "s * m^-1" => u"s * m^-1",
    "dia" => u"d",
    "kilograms" => u"kg",
    "kg * s^-1" => u"kg * s^-1",
    "kg * m^-3" => u"kg * m^-3",
    "kg * kmol^-1"=> u"kg * kmol^-1",
    "kg_water * kg_air^-1" => u"kg * kg^-1",
    "kg_vapour * J^-1" => u"kg * J^-1",
    "g * mol_CH2O^-1" => u"g * mol^-1",
    "mg_CO2 * J^-1" => u"mg * J^-1",
    "mg * s^-1" => u"mg * s^-1",
    "J * kg_water^-1" => u"J * kg^-1",
    "J * kmol^-1 * K^-1" => u"J * kmol^-1 * K^-1",
    "J * m^-3 * K^-1" => u"J * m^-3 * K^-1",
    "J * K^-1 * m^-2" => u"J * K^-1 * m^-2",
    "J * K^-1 * kg^-1" => u"J * K^-1 * kg^-1",
    "MXN * kg" => u"MXN * kg",
    "EUR * kW^-1 * h^-1" => u"EUR * kW^-1 * h^-1",
    "EUR * kg^-1" => u"EUR * kg^-1",
    "EURos / m^2" => u"EUR * m^-2",
    "MXN * m^-2" => u"MXN * m^-2",
    "MXN * m^-2 * s^-1" => u"MXN * m^-2 * s^-1",
    "MXN *m^-3" => u"MXN *m^-3",
    "W" => u"W",
    "kW * h * m^-3" => u"kW * h * m^-3",
    "W * m^-1 * K^-1" => u"W * m^-1 * K^-1",
    "W * m^-2" => u"W * m^-2",
    "W * m^-2 * K^-1" => u"W * m^-2 * K^-1",
    "W * m^-2 * K^-4" => u"W * m^-2 * K^-4",
    "W * m_cover^-2 * K^-1" => u"W * m^-2 * K^-1",
    "Pa^-2" => u"Pa^-2",
    "Pa * K^-1" => u"Pa * K^-1",
    "ppm^-2" => u"ppm^-2",
    "ppm * mg^-1 * m^3" => u"ppm * mg^-1 * m^3",
    "m^3 * m^-2 * s^-1" => u"m^3 * m^-2 * s^-1",

    # Add more unit mappings as needed based on the content of your JSON file

    "m^3 * g^-1 * s^-1 * mol_CO2 * mol_air^-1" => u"m^3 * g^-1 * s^-1 * mol * mol^-1",
    "mol_phot * m^-2 * d^-1" => u"mol * m^-2 * d^-1",
    "mol_phot * m^-2 * s^-1" => u"mol * m^-2 * s^-1",
    "mol * m^-2 * s^-1" => u"mol * m^-2 * s^-1",
    "mol_CO2 * mol_air^-1" => u"mol * mol^-1",
    "mol_air * mol_CO2^-1" => u"mol * mol^-1",
    "m^2 * d * mol_phot^-1" => u"m^2 * d * mol^-1",
    "Pa^-1" => u"Pa^-1",
    "s * m^-1" => u"s * m^-1",
    "g * d * mol_CO2^-1" => u"g * d * mol^-1",
    "g * s * mol_CO2^-1" => u"g * s * mol^-1",
    "mol_O2 * mol_CO2^-1" => u"mol * mol^-1",
    "mol_O2 * mol_air^-1" => u"mol * mol^-1",
    "mol_CO2 * m^-2 * d^-1" => u"mol * m^-2 * d^-1",
    "mol_CO2 * m^-2 * d^-1" => u"mol * m^-2 * s^-1",
    "mol_O2 * mol_air^-1" => u"mol * mol^-1",
    "mol_CO2 * m^-2 * s^-1" => u"mol * m^-2 * s^-1",
    "mol * mol^-1" => u"mol * mol^-1",
    "Pa * Pa^-1" => u"Pa * Pa^-1",
    "m^2 * kg * s^-2 * K^-1 * mol^-1" => u"m^2 * kg * s^-2 * K^-1 * mol^-1",
    "kg * m^2 * s^-2 * mol^-1" => u"kg * m^2 * s^-2 * mol^-1",
    "m^2 * s * mol_phot^-1" => u"m^2 * s * mol^-1"
    )

"""
    importdictf_nounits(const_dict)

Function to import constants from a dictionary and create global variables
# Arguments:
- const_dict::Dict{String, Dict}: A dictionary containing constants and their metadata. 
   Each key represents the variable name, and the associated value is a dictionary with the following fields:
     - "val": The numeric value of the constant.
     - "units": The string representation of the unit.
     - "desc": A brief description of the constant.

# Requirements:
 - The `unit_mapping` global variable should be populated with valid unit mappings.

# Behavior:
 - For each key-value pair in the dictionary:
     - If the required fields ("val", "units", "desc") exist:
         - Maps the unit string to a `Unitful.jl` unit using `unit_mapping`.
         - Dynamically creates global variables for:
             - The constant's value (without units).
             - The unit of the constant (if valid in `unit_mapping`).
             - A description of the constant.
     - Logs an error if:
         - Any of the required fields are missing.
         - The unit is not found in `unit_mapping`.
         - An unexpected error occurs during processing.
 - If the dictionary is empty, it logs a warning message.

# Notes:
 - This function uses `eval` to create global variables dynamically. Use with caution.
 - Global variables are named after the keys in the dictionary, with `_units` and `_description` appended for units and descriptions respectively.
"""
function importdictf_nounits(const_dict)
    if isempty(const_dict)
        println("The dictionary is empty or couldn't be read correctly.")
    else
        # Iterate over the dictionary and create global variables with units
        for (key, value) in const_dict
            try
                # Check if the key contains the expected structure
                if haskey(value, "val") && haskey(value, "units") && haskey(value, "desc")
                    var_value = value["val"]
                    var_units = value["units"]
                    var_description = value["desc"]

                    # Map the unit string from JSON to a Unitful unit
                    if haskey(unit_mapping, var_units)
                        unitful_unit = unit_mapping[var_units]
                        unitless_value = var_value 

                        # Dynamically create a global variable with the name from the JSON key
                        global_var_name = Symbol(key)
                        eval(:(global $global_var_name = $unitless_value))  # Ensure the unitful_value is evaluated globally

                        # Optionally, create global variables for units and description
                        units_var_name = Symbol(key * "_units")
                        eval(:(global $units_var_name = $unitful_unit))

                        description_var_name = Symbol(key * "_description")
                        eval(:(global $description_var_name = $var_description))

                    else
                        println("Error: Unit '$var_units' not found in unit mapping for key $key.")
                    end
                else
                    println("Error: Missing 'val', 'units', or 'desc' field for key $key.")
                end

            catch e
                println("Error processing key $key: ", e)
            end
        end
    end
end


# Define a global structure to hold constants for better organization
mutable struct Constants
    climate::Dict
    crop_photo::Dict
end

# Function to read and parse a JSON file, returning a dictionary
function read_json_file(file_path::String)::Dict
    try
        return JSON.parsefile(file_path)
    catch e
        println("Error reading JSON file: ", e)
        return Dict()  # Return an empty dictionary on error
    end
end

# Function to process the dictionary and import constants into the program
function import_constants(const_dict::Dict)
    importdictf_nounits(const_dict)
end

# Main function to load and import constants from multiple files
function load_constants()
    # Define paths for climate and crop photo constants
    
    climate_file_path = "src/Models/parameters/climate_parameters.json"
    crop_photo_file_path = "src/Models/parameters/photosynthesis_parameters.json"

    # Load climate constants from the JSON file
    climate_constants = read_json_file(climate_file_path)
    
    # Load crop photo constants from the second JSON file
    crop_photo_constants = read_json_file(crop_photo_file_path)
    
    # Create an instance of Constants with the loaded dictionaries
    constants = Constants(climate_constants, crop_photo_constants)
    
    # Import constants into the program (process them without units)
    import_constants(constants.climate)
    import_constants(constants.crop_photo)
    
    return constants
end


# Call the main function to load and process constants
constants = load_constants()


