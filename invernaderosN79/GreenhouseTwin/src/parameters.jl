# =============================================================================
# parameters.jl -- one flat, concrete NamedTuple of Float64 parameters per model part
# =============================================================================
# Parameter files (configfiles/): climate_n79.json (climate + photosynthesis + noise), crop_n105.json (crop + noise).
# Each entry: {"val", "units", "group", "status" (fixed | inferred), "prior" ([min, centre, max], inferred only),
# "desc", "source"}. For inferred parameters "val" is the prior centre.

"""
    load_params(json_path) -> NamedTuple

All parameter values of one file as a concrete NamedTuple of Float64 (fast field access inside the ODE).
"""
function load_params(json_path::AbstractString)
    raw = JSON.parsefile(json_path)
    d = Dict{Symbol,Float64}()
    for (k, v) in raw
        d[Symbol(k)] = Float64(v["val"])
    end
    return (; sort(collect(d); by = first)...)
end

"Raw metadata Dict of a parameter file (units, group, status, prior, desc, source)."
param_meta(json_path::AbstractString) = JSON.parsefile(json_path)

"Names of the parameters with status \"inferred\" (sorted)."
inferred_names(json_path::AbstractString) =
    sort([k for (k, v) in param_meta(json_path) if get(v, "status", "") == "inferred"])

"Dict name => [min, centre, max] of the inferred parameters."
priors_of(json_path::AbstractString) =
    Dict(k => Float64.(v["prior"]) for (k, v) in param_meta(json_path) if haskey(v, "prior"))

"""
    update_params(base::NamedTuple, names, x) -> NamedTuple

Copy of `base` with the parameters in `names` set to `x` (same order). Every name must already exist in `base`
(a misspelt name is an error, not a silently added field).
"""
function update_params(base::NamedTuple, names, x)
    syms = Tuple(Symbol(n) for n in names)
    for s in syms
        hasproperty(base, s) || error("update_params: unknown parameter $s")
    end
    vals = Tuple(Float64(v) for v in x)
    return merge(base, NamedTuple{syms}(vals))
end
