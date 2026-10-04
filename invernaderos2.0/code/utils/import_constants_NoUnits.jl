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

    
# Ensure const_dict is correctly populated before proceeding
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



# Fix the loading blocks
base_dir = pwd()
json_file_path = joinpath(base_dir, "configfiles", "constants_climate.json")

try
    global const_dict = JSON.parsefile(json_file_path)
    importdictf_nounits(const_dict)
catch e
    println("Error reading JSON file $json_file_path: ", e)
end
 """

# Build paths using joinpath (this works on Mac, Linux, and Windows automatically)
file_path = joinpath(base_dir, "configfiles", "constants_climate.json")

# Define the path to the JSON file
#json_file_path = "/Users/capella/Documents/trabajo/academia/unam/proyectos investigacion/bio/invernaderos/invernaderos2.0/code/configfiles/constants_climate.json"

# Read and parse the JSON file
const_dict = Dict()  # Initialize as an empty dictionary to avoid referencing issues
try
    global const_dict = JSON.parsefile(json_file_path)  # Parse the JSON file into a dictionary
catch e
    println("Error reading JSON file: ", e)
end
importdictf_nounits(const_dict)
"""
# Example usage:



# Build paths using joinpath (this works on Mac, Linux, and Windows automatically)

# Fix the loading blocks
base_dir = pwd()
json_file_path = joinpath(base_dir, "configfiles", "constants_cropphoto.json")

try
    global const_dict = JSON.parsefile(json_file_path)
    importdictf_nounits(const_dict)
catch e
    println("Error reading JSON file $json_file_path: ", e)
end

"""
file_path = joinpath(base_dir, "configfiles", "constants_cropphoto.json")

# Read and parse the JSON file
const_dict = Dict()  # Initialize as an empty dictionary to avoid referencing issues
try
    global const_dict = JSON.parsefile(json_file_path)  # Parse the JSON file into a dictionary
catch e
    println("Error reading JSON file: ", e)
end
importdictf_nounits(const_dict)
"""