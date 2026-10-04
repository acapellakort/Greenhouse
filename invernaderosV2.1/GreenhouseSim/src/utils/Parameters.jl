# src/utils/Parameters.jl

# --- 1. Unit Definitions ---
# We check if we are in a module to avoid re-registering during precompilation
Unitful.register(@__MODULE__)

# Use 'const' for units to prevent redefinition errors
@unit USD "USD" Dollar 1.0 false
@unit MXN "MXN" Pesos 1.0 false
@unit EUR "EUR" Euros 1.0 false
# 'ppm' is often already in Unitful; if you get a warning, you can skip this one
# @unit ppm "ppm" part_per_million 1e-6 false 

const unit_mapping = Dict{String, Any}(
    "1" => u"NoUnits",
    "m" => u"m",
    "meters" => u"m",
    "metros" => u"m",
    "J * K^-1 * m^-2" => u"J * K^-1 * m^-2",
    "W * m^-2 * K^-1" => u"W * m^-2 * K^-1",
    "C" => u"°C",
    # ... (Keep your other mapping entries here) ...
    "m^2 * s * mol_phot^-1" => u"m^2 * s * mol^-1"
)

# --- 2. Struct Definitions ---
struct ModelParams
    vals::NamedTuple
    units::Dict{Symbol, Any}
    meta::Dict{Symbol, String}
end

# Property access: p.alpha1
Base.getproperty(p::ModelParams, s::Symbol) = s in fieldnames(ModelParams) ? getfield(p, s) : p.vals[s]

# --- 3. Functions ---
function load_model_params(json_path::String, mapping::Dict = unit_mapping)
    raw = JSON.parsefile(json_path)
    vals_dict = Dict{Symbol, Float64}()
    meta_dict = Dict{Symbol, String}()
    units_dict = Dict{Symbol, Any}()
    
    for (k, v) in raw
        sym = Symbol(k)
        vals_dict[sym] = Float64(v["val"])
        meta_dict[sym] = v["desc"]
        u_str = v["units"]
        units_dict[sym] = get(mapping, u_str, u"NoUnits")
    end
    
    return ModelParams(NamedTuple(vals_dict), units_dict, meta_dict)
end

function merge_params(p1::ModelParams, p2::ModelParams)
    new_vals = merge(p1.vals, p2.vals)
    new_units = merge(p1.units, p2.units)
    new_meta = merge(p1.meta, p2.meta)
    return ModelParams(new_vals, new_units, new_meta)
end