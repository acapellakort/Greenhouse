# =============================================================================
# parameters.jl  --  concrete, type-stable parameter container
# =============================================================================
#
# All model constants live in ONE flat NamedTuple of Float64 values. A
# NamedTuple with homogeneous concrete element types is a fully concrete type,
# so `p.alpha1` / `(; alpha1) = p` is inferable and zero-cost inside the RHS.
#
# This REPLACES 2.0/V2.1's `import_constants_NoUnits.jl`, which defined every
# constant as a non-`const` global via `eval(:(global ...))` -- the single
# biggest source of type instability and the ~13 GiB of allocations reported
# in the old README. There are no globals here.

"""
    load_params(climate_json, crop_json) -> NamedTuple

Read the two JSON config files (same format as 2.0: each entry has "val",
"units", "desc") and return a single flat `NamedTuple` mapping every constant
name to its `Float64` value. Climate and crop constants are merged; if a name
appears in both, the crop file wins (none currently overlap).
"""
function load_params(climate_json::AbstractString, crop_json::AbstractString)
    raw = merge(JSON.parsefile(climate_json), JSON.parsefile(crop_json))
    d = Dict{Symbol,Float64}()
    for (k, v) in raw
        # Every entry stores its numeric value under "val"; coerce to Float64
        # so the resulting NamedTuple is homogeneous (fastest concrete case).
        d[Symbol(k)] = Float64(v["val"])
    end
    return (; d...)   # Dict{Symbol,Float64}  ->  concrete NamedTuple
end

"""
    update_params(base::NamedTuple, names, x) -> NamedTuple

Return a copy of `base` with the parameters in `names` set to the values in
`x` (same order). `names` may be `Symbol`s or `String`s and is chosen at run
time -- this is how you select which parameters to infer, replacing 2.0's
global-variable trick with no loss of flexibility.

    p = update_params(base, ["beta2", "gamma3"], [0.7, 275.0])

The returned NamedTuple has exactly the same fields as `base` (so it stays a
concrete type); only the listed values change. Call this ONCE per MCMC
proposal, then pass the result into `solve` (a function barrier), so the ODE
inner loop always sees a concretely-typed parameter set.
"""
function update_params(base::NamedTuple, names, x)
    syms = Tuple(Symbol(n) for n in names)
    vals = Tuple(Float64(v) for v in x)
    return merge(base, NamedTuple{syms}(vals))
end
