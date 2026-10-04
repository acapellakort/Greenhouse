# =============================================================================
# measured.jl  --  drive growth from MEASURED greenhouse climate (for calibration)
# =============================================================================
#
# For calibrating the growth parameters we drive photosynthesis directly from the
# measured greenhouse climate (canopy temperature, CO2, light) instead of
# re-simulating the climate ODE. This isolates the growth parameters from any
# climate-model error -- standard practice for crop-model calibration -- and is
# very fast (no ODE; a whole season is ~30k assimilation evaluations).
#
# Data: Autonomous Greenhouse Challenge 2018 (cucumber). Per team:
#   Greenhouse_climate.csv : GHtime (Excel serial), Tair [C], CO2air [ppm], AssimLight [%]
#   meteo.csv              : Iglob [W/m2] (outside), row-aligned with the climate file
#   Production.csv         : time (Excel serial), Total_Prod_cum [kg FW/m2]

"Season of measured climate, grouped by day. Each day: 3xN matrix [Tair_K; CO2_mg/m3; Iw]."
struct MeasuredSeason
    ndays::Int
    day_samples::Vector{Matrix{Float64}}
    dt::Float64                    # seconds between samples (300 = 5 min)
end

# Excel serial date -> day index relative to a season start serial.
_excel_day(serial, day0) = Int(floor(serial - day0))

# ppm CO2 -> mg/m3 at temperature Tc [C] (ideal gas): mg/m3 = ppm * 536.4 / T[K].
_co2_ppm_to_mgm3(ppm, Tc) = ppm * 536.4 / (Tc + 273.15)

"""
    load_measured_climate(gh_csv, meteo_csv, p; day0_serial, dt) -> MeasuredSeason

Build a per-day sample season from the measured greenhouse climate. `p` is the
climate+crop parameter NamedTuple; its radiation-transmission constants
(`eta1, tau1, eta2, alpha12`) map outside `Iglob` + lamp `%` to canopy light `Iw`
using the same convention as the model RHS.
"""
function load_measured_climate(gh_csv::AbstractString, meteo_csv::AbstractString, p::NamedTuple;
                               day0_serial::Float64 = 43326.0, dt::Float64 = 300.0)
    (; eta1, tau1, eta2, alpha12) = p
    trans = (1 - eta1) * tau1 * eta2

    gh = DataFrame(CSV.File(gh_csv))
    mt = DataFrame(CSV.File(meteo_csv))
    n  = min(nrow(gh), nrow(mt))

    buckets = Dict{Int, Vector{NTuple{3,Float64}}}()
    isnum(x) = (x isa Number) && !isnan(x)

    for i in 1:n
        g = gh.GHtime[i]
        isnum(g) || continue
        Tc = gh.Tair[i]; Co = gh.CO2air[i]; La = gh.AssimLight[i]; Ig = mt.Iglob[i]
        (isnum(Tc) && isnum(Co)) || continue          # need canopy T and CO2
        Igv = isnum(Ig) ? Float64(Ig) : 0.0
        Lav = isnum(La) ? Float64(La) : 0.0
        d = _excel_day(g, day0_serial)
        d >= 0 || continue
        Tk  = Tc + 273.15
        Cmg = _co2_ppm_to_mgm3(Co, Tc)
        Iw  = trans * Igv + alpha12 * (Lav / 100.0)    # canopy light (model units)
        push!(get!(buckets, d, NTuple{3,Float64}[]), (Tk, Cmg, Iw))
    end

    ndays = isempty(buckets) ? 0 : maximum(keys(buckets)) + 1
    day_samples = Vector{Matrix{Float64}}(undef, ndays)
    for d in 0:(ndays - 1)
        v = get(buckets, d, NTuple{3,Float64}[])
        M = Matrix{Float64}(undef, 3, length(v))
        for (j, s) in enumerate(v)
            M[1, j] = s[1]; M[2, j] = s[2]; M[3, j] = s[3]
        end
        day_samples[d + 1] = M
    end
    return MeasuredSeason(ndays, day_samples, dt)
end

"Daily gross assimilate [g CH2O m^-2 d^-1] from one day's measured samples at `lai`."
function daily_assimilation_measured(M::Matrix{Float64}, lai::Float64, p::NamedTuple, dt::Float64,
                                     k_ext::Float64 = K_EXT)
    size(M, 2) == 0 && return 0.0
    s = 0.0
    @inbounds for j in 1:size(M, 2)
        s += assimilation(M[1, j], M[2, j], M[3, j], lai_eff(lai, k_ext), p)   # canopy interception
    end
    return s * dt * 1e-3 * (30.0 / 44.0)                        # -> g CH2O m^-2 d^-1
end

"Mean canopy temperature [deg C] for a day's samples (fallback 20 C if empty)."
function daily_mean_T(M::Matrix{Float64})
    size(M, 2) == 0 && return 20.0
    return sum(@view M[1, :]) / size(M, 2) - 273.15
end

"""
    simulate_growth_measured(p, gp, season; LAI0) -> NamedTuple of daily series

Season simulation of the cucumber growth model driven by measured climate.
Returns day, LAI, yield_FW [kg/m2], W_fruit, Pg, and the final state.
"""
function simulate_growth_measured(p::NamedTuple, gp, season::MeasuredSeason; LAI0 = 0.5)
    s = init_growth_state(gp; LAI0 = LAI0)
    day = Int[]; LAI = Float64[]; yieldFW = Float64[]; Wfruit = Float64[]; Pg = Float64[]; node = Float64[]
    for d in 0:(season.ndays - 1)
        M = season.day_samples[d + 1]
        Pgd = daily_assimilation_measured(M, s.LAI, p, season.dt, gp.k_ext)
        grow!(s, Pgd, daily_mean_T(M), gp, d)
        push!(day, d); push!(LAI, s.LAI); push!(yieldFW, s.yield_FW)
        push!(Wfruit, isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits))
        push!(Pg, Pgd); push!(node, s.node)
    end
    return (; day, LAI, yield_FW = yieldFW, W_fruit = Wfruit, Pg, node, state = s)
end

"""
    load_cropmanagement(csv_path; leaf_area, stem_density, N_window, day0_date)
        -> (days, N_leaves, LAI_obs)

Read CropManagement.csv (`weeks`, `N_leaves` = cumulative leaves per stem) and
build two observables aligned to the season day index:
  * `N_leaves`  — cumulative leaf (≈ node) number per stem, for development;
  * `LAI_obs`   — a DERIVED standing LAI = min(N_leaves, N_window)·leaf_area·stem_density.

The three conversion constants are literature assumptions (no measured LAI /
plant density in the dataset) — adjust `leaf_area`, `stem_density`, `N_window`
if the real values are known.
"""
function load_cropmanagement(csv_path::AbstractString;
                             leaf_area::Float64    = 0.05,   # m^2 per leaf
                             stem_density::Float64 = 2.5,    # stems / m^2  (ASSUMED)
                             N_window::Float64     = 20.0,   # functional leaves kept per stem
                             day0_date::Date       = Date(2018, 8, 14))
    df = DataFrame(CSV.File(csv_path))
    days    = [Dates.value((Date(2018, 1, 1) + Day((Int(w) - 1) * 7)) - day0_date) for w in df.weeks]
    Nleaves = Float64.(df.N_leaves)
    keep    = .!isnan.(Nleaves) .& (days .>= 0)
    d  = days[keep]; Nl = Nleaves[keep]
    standing = min.(Nl, N_window)
    LAI_obs  = standing .* leaf_area .* stem_density
    return d, Nl, LAI_obs
end

"""
    load_production(prod_csv; day0_serial) -> (days, cum_yield)

Observed cumulative cucumber yield [kg FW/m2] and its day index (day 0 = day0_serial).
"""
function load_production(prod_csv::AbstractString; day0_serial::Float64 = 43326.0)
    df = DataFrame(CSV.File(prod_csv))
    days = _excel_day.(Float64.(df.time), day0_serial)
    return days, Float64.(df.Total_Prod_cum)
end
