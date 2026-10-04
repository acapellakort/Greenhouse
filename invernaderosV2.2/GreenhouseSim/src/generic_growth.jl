# =============================================================================
# generic_growth.jl  --  Marcelis/Heuvelink generic dry-matter partitioning
# =============================================================================
#
# Julia port of Aaron Velez's Generic_Gh_Model (Python) which unifies:
#   M94: Marcelis (1994), Annals of Botany 74, 43-52.   (cucumber)
#   H96: Heuvelink (1996), Annals of Botany 77, 71-80.  (tomato / TOMSIM)
#
# This file is DEFAULT-OFF: grow!() and grow_cohort!() pass
# `generic_params = nothing` by default, preserving every existing
# V2.2 result exactly. Pass a GenericParameters instance to activate
# the physics here.
#
# LAI, photosynthesis, and CO2 are external to this model (unchanged in V2.2).
# The generic model covers:
#   - reproductive event appearance (flowers / trusses per plant per day)
#   - organ development (normalized thermal clock 0→1)
#   - organ sink strength (g DM organ^-1 d^-1)
#   - vegetative demand (g DM plant^-1 d^-1)
#   - daily assimilate allocation (M94 eq.4 / H96 eq.1)
# =============================================================================

# ---------------------------------------------------------------------------
# Parameter container
# ---------------------------------------------------------------------------

"""
    GenericParameters

Immutable struct holding all coefficients for one crop parameterization.
Covers reproductive appearance, thermal development clock, organ sink,
and vegetative demand. All fields are Float64. Domain tuples are (lo, hi)
in °C or [0,1].

Construct with keyword arguments (Base.@kwdef):

    gp = GenericParameters(s=4.675, p=91.9009, peak=0.4764, ...)

See `CUCUMBER_PARAMS` and `TOMATO_PARAMS` for calibrated instances.
"""
Base.@kwdef struct GenericParameters
    # Sink-strength Bell curve (generic_sink)
    s    :: Float64    # shape × endpoint  (dimensionless)
    p    :: Float64    # exponent (dimensionless)
    peak :: Float64    # normalized stage at peak [0,1]
    A0   :: Float64    # amplitude intercept [g DM d^-1]
    AT   :: Float64    # amplitude temperature slope [g DM °C^-1 d^-1]
    # Thermal development clock (generic_development_rate)
    kappa :: Float64   # rate constant [stage °C^-alpha d^-1]
    alpha :: Float64   # temperature exponent (1 = linear, <1 = sublinear)
    base  :: Float64   # base temperature [°C]
    # Reproductive appearance rate (generic_appearance)
    beta0 :: Float64   # intercept [events plant^-1 d^-1]
    beta1 :: Float64   # slope [events plant^-1 °C^-1 d^-1]
    # Vegetative demand (generic_vegetation)
    V20   :: Float64   # demand at 20 °C [g DM plant^-1 d^-1]
    VT    :: Float64   # temperature slope [g DM plant^-1 °C^-1 d^-1]
    # Valid domains (lo, hi) in °C, or [0,1] for stage
    clock_domain       :: Tuple{Float64,Float64}
    appearance_domain  :: Tuple{Float64,Float64}
    sink_domain        :: Tuple{Float64,Float64}
    vegetation_domain  :: Tuple{Float64,Float64}
end

# ---------------------------------------------------------------------------
# Calibrated instances (exact coefficients from models.py)
# ---------------------------------------------------------------------------

"""
    CUCUMBER_PARAMS :: GenericParameters

Cucumber parameterization from Marcelis (1994) cast into the generic form.
Appearance rate fitted at PAR = 6.8 MJ m^-2 d^-1.
Thermal endpoint: 275 degree-C days above 10 °C.
"""
const _PAR_CUCUMBER = 6.8
const _CUCUMBER_ENDPOINT_DD = 275.0
const _par_factor = 1.0 - exp(-0.5 - 0.5 * _PAR_CUCUMBER)

const CUCUMBER_PARAMS = GenericParameters(
    s    = 0.0170 * _CUCUMBER_ENDPOINT_DD,                 # 4.675
    p    = 1.0 + 1.0 / 0.0111,                            # 91.1261...
    peak = 131.0 / _CUCUMBER_ENDPOINT_DD,                  # 0.47636...
    A0   = -10.0 * 0.0170 * 60.7 / 0.0111,                # -927.027...
    AT   =         0.0170 * 60.7 / 0.0111,                 # 92.7027...
    kappa = 1.0 / _CUCUMBER_ENDPOINT_DD,                   # 0.003636...
    alpha = 1.0,
    base  = 10.0,
    beta0 = -0.75 * _par_factor,
    beta1 =  0.09 * _par_factor,
    V20   = 8.3,
    VT    = 0.25,
    clock_domain      = (10.0, 30.0),
    appearance_domain = (17.5, 30.0),
    sink_domain       = (17.5, 30.0),
    vegetation_domain = (18.0, 24.0),
)

"""
    TOMATO_PARAMS :: GenericParameters

Tomato parameterization from Heuvelink (1996) / TOMSIM cast into the generic form.
Thermal clock fitted numerically to the TOMSIM recurrence.
"""
const TOMATO_PARAMS = GenericParameters(
    s    = 4.34,
    p    = 1.31 / (1.31 - 1.0),                            # 4.241379...
    peak = 0.278 - log(1.31 - 1.0) / 4.34,                 # 0.561416...
    A0   = 0.138 * 4.34 / (1.31 - 1.0),                    # 0.585103...
    AT   = 0.0,
    kappa = 0.0072996177216267171,
    alpha = 0.45412065794525647,
    base  = 12.822328748552815,
    beta0 = 0.020768970820626603,
    beta1 = 0.0062165461513216625,
    V20   = 2.8,
    VT    = 0.0,
    clock_domain      = (17.5, 30.0),
    appearance_domain = (17.5, 30.0),
    sink_domain       = (18.0, 24.0),
    vegetation_domain = (18.0, 24.0),
)

# ---------------------------------------------------------------------------
# Generic functions (unified, work for both crops)
# ---------------------------------------------------------------------------

"""
    generic_appearance(T, gp) -> events plant^-1 d^-1

Affine reproductive appearance rate (flowers for cucumber, trusses for tomato).
Returns 0 below the appearance domain and at the lower boundary.
"""
function generic_appearance(T::Float64, gp::GenericParameters)
    lo, hi = gp.appearance_domain
    T < lo && return 0.0
    T > hi && (T = hi)
    return max(gp.beta0 + gp.beta1 * T, 0.0)
end

"""
    generic_development_rate(T, gp) -> stage d^-1

Normalized thermal clock: kappa * max(T - base, 0)^alpha.
Returns 0 below the base temperature and outside the clock domain.
"""
function generic_development_rate(T::Float64, gp::GenericParameters)
    lo, hi = gp.clock_domain
    T = clamp(T, lo, hi)
    excess = T - gp.base
    excess <= 0.0 && return 0.0
    return max(gp.kappa * excess ^ gp.alpha, 0.0)
end

"""
    generic_sink(T, z, n, gp) -> g DM d^-1

Generic Bell-curve organ sink strength for a cohort of `n` organs at
normalized development stage `z ∈ [0,1]` and temperature `T`.
`n` may be Float64 (fruits/m²) — the formula is linear in `n`.

Returns 0 outside the sink domain or for z outside [0,1].
"""
function generic_sink(T::Float64, z::Float64, n::Real, gp::GenericParameters)
    lo, hi = gp.sink_domain
    (T < lo || T > hi || z < 0.0 || z > 1.0) && return 0.0
    amplitude = max(gp.A0 + gp.AT * T, 0.0)
    amplitude == 0.0 && return 0.0
    u = exp(-gp.s * (z - gp.peak)) / (gp.p - 1.0)
    return max(Float64(n) * amplitude * u / (1.0 + u)^gp.p, 0.0)
end

"""
    generic_vegetation(T, gp) -> g DM plant^-1 d^-1

Aggregate vegetative demand (V20 + VT*(T-20)).
Clamped to the vegetation domain and to non-negative values.
"""
function generic_vegetation(T::Float64, gp::GenericParameters)
    lo, hi = gp.vegetation_domain
    T = clamp(T, lo, hi)
    return max(gp.V20 + gp.VT * (T - 20.0), 0.0)
end

"""
    allocate_daily(potentials, supply, reserve=0.0) -> (growths, next_reserve)

M94 eq.4 / H96 eq.1 with daily reserve carry-over.
`potentials`: iterable of non-negative demand amounts [g DM].
`supply`:     daily assimilate available [g DM].
`reserve`:    carry-over from previous day [g DM].
Returns a `Vector{Float64}` of per-sink growth amounts and the leftover reserve.

When total demand is zero, all available mass enters the reserve.
"""
function allocate_daily(potentials, supply::Float64, reserve::Float64 = 0.0)
    demands   = collect(Float64, potentials)
    available = max(supply, 0.0) + max(reserve, 0.0)
    total     = sum(demands)
    if total <= 0.0
        return zeros(length(demands)), available   # nothing to grow; bank it all
    end
    growth = demands .* (min(total, available) / total)
    remaining = available - sum(growth)
    remaining < 0.0 && (remaining = 0.0)
    return growth, remaining
end

# ---------------------------------------------------------------------------
# Original crop-specific reference functions (M94 / H96 formulas verbatim)
# Kept for validation: generic_* should reproduce these within floating-point.
# ---------------------------------------------------------------------------

"""
    cucumber_appearance_m94(T; PAR=6.8) -> flowers plant^-1 d^-1

M94 eq.7, p.46. Fitted at PAR = 6.8 MJ m^-2 d^-1.
Domain: T ∈ [17.5, 30.0] °C.
"""
function cucumber_appearance_m94(T::Float64; PAR::Float64 = 6.8)
    T = clamp(T, 17.5, 30.0)
    return max((-0.75 + 0.09 * T) * (1.0 - exp(-0.5 - 0.5 * PAR)), 0.0)
end

"""
    cucumber_sink_m94(T, thermal_age) -> g DM fruit^-1 d^-1

M94 eq.6, p.45. Bell sigmoid in thermal age [°C·d] above 10 °C.
Domain: T ∈ [17.5, 30.0] °C, thermal_age ∈ [0, 270] °C·d.
"""
function cucumber_sink_m94(T::Float64, thermal_age::Float64)
    T = clamp(T, 17.5, 30.0)
    x = clamp(thermal_age, 0.0, 270.0)
    e = exp(-0.0170 * (x - 131.0))
    return (T - 10.0) * 0.0170 * 60.7 * e / (1.0 + 0.0111 * e)^(1.0 + 1.0/0.0111)
end

"""
    cucumber_vegetation_m94(T) -> g DM plant^-1 d^-1

M94 p.45 unnumbered linear rule; includes roots.
Domain: T ∈ [18.0, 24.0] °C.
"""
function cucumber_vegetation_m94(T::Float64)
    T = clamp(T, 18.0, 24.0)
    return 7.8 + (9.3 - 7.8) * (T - 18.0) / 6.0
end

"""
    tomato_appearance_h96(T) -> trusses plant^-1 d^-1

H96 eq.2, p.72.
Domain: T > 0 °C.
"""
function tomato_appearance_h96(T::Float64)
    T <= 0.0 && return 0.0
    return max(-0.2903 + 0.1454 * log(T), 0.0)
end

"""
    tomato_development_rate_h96(T, stage) -> stage d^-1

H96 eq.3, p.73. Stage is the PREVIOUS day's normalized development.
Domain: T > 0 °C, stage ∈ [0,1].
"""
function tomato_development_rate_h96(T::Float64, stage::Float64)
    T <= 0.0 && return 0.0
    z = clamp(stage, 0.0, 1.0)
    poly = 0.0392 - 0.213*z + 0.451*z^2 - 0.240*z^3
    return max(0.0181 + log(T / 20.0) * poly, 0.0)
end

"""
    tomato_sink_h96(stage, fruits) -> g DM truss^-1 d^-1

H96 eq.5, p.73. `fruits` is the number of fruits per truss (integer ≥ 1).
"""
function tomato_sink_h96(stage::Float64, fruits::Int)
    fruits < 1 && error("fruits must be ≥ 1")
    z = clamp(stage, 0.0, 1.0)
    a, b, c, d = 0.138, 4.34, 0.278, 1.31
    num = a * b * (1.0 + exp(-b * (z - c)))^(1.0 / (1.0 - d))
    den = (d - 1.0) * (exp(b * (z - c)) + 1.0)
    return max(fruits * num / den, 0.0)
end

"""
    tomato_vegetation_h96() -> g DM plant^-1 d^-1

H96 p.73 unnumbered constant: above-ground vegetative demand.
"""
tomato_vegetation_h96() = 2.8
