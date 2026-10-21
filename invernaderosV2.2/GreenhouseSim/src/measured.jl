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
function simulate_growth_measured(p::NamedTuple, gp, season::MeasuredSeason; LAI0 = 0.5,
                                  generic_params::Union{Nothing, GenericParameters} = nothing,
                                  fruit_keep = nothing)
    # Note 69: fruit_keep = per-day vector (index d+1) of the grower's kept-fruit share (generic path); harvested fruit
    # numbers and their mean age (days since set) are tracked per day (harv_n, harv_age) for the FruitGrw observation.
    s = init_growth_state(gp; LAI0 = LAI0)
    if hasproperty(gp, :gen_canopy) && gp.gen_canopy > 0.5      # Note 94: transplant state for the canopy path
        s.node = gp.gen_N0; s.leaves = floor(gp.gen_N0)
        s.W_leaf = LAI0 / seasonal_sla_m2g(gp.gen_doy0, gp)    # leaf DW behind the supplied initial LAI
        nl = s.W_leaf * (1.0 - gp.gen_fL) / gp.gen_fL           # non-leaf DW in the vegetative split
        fst = gp.frac_stem / (gp.frac_stem + gp.frac_root); s.W_stem = fst * nl; s.W_root = (1.0 - fst) * nl
    end
    day = Int[]; LAI = Float64[]; yieldFW = Float64[]; Wfruit = Float64[]; Pg = Float64[]; node = Float64[]
    harv_n = Float64[]; harv_age = Float64[]
    Wl = Float64[]; Ws = Float64[]; Wr = Float64[]; Ydw = Float64[]     # Note 79: organ masses for respiration
    Lv = Float64[]                                                       # Note 94: retained leaves per stem
    live = Dict{UInt,Tuple{Int,Fruit}}()
    for d in 0:(season.ndays - 1)
        M = season.day_samples[d + 1]
        Pgd = daily_assimilation_measured(M, s.LAI, p, season.dt, gp.k_ext)
        parMJ = size(M, 2) == 0 ? NaN : sum(@view M[3, :]) * season.dt * 1e-6      # daily canopy PAR [MJ m-2 d-1]
        keep = (fruit_keep === nothing) ? NaN : Float64(fruit_keep[min(d + 1, length(fruit_keep))])
        grow!(s, Pgd, daily_mean_T(M), gp, d; generic_params = generic_params, fruit_keep = keep, par_MJ = parMJ)
        cur = Set(objectid(f) for f in s.fruits); hn = 0.0; ha = 0.0
        for (k, (sd, f)) in collect(live)
            if !(k in cur); hn += f.n; ha += f.n * (d - sd); delete!(live, k); end
        end
        for f in s.fruits
            k = objectid(f); (!haskey(live, k) && f.n > 0.0) && (live[k] = (d, f))
        end
        push!(day, d); push!(LAI, s.LAI); push!(yieldFW, s.yield_FW)
        push!(Wfruit, isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits))
        push!(Pg, Pgd); push!(node, s.node); push!(harv_n, hn); push!(harv_age, hn > 0 ? ha / hn : NaN)
        push!(Wl, s.W_leaf); push!(Ws, s.W_stem); push!(Wr, s.W_root); push!(Ydw, s.yield_DW); push!(Lv, s.leaves)
    end
    return (; day, LAI, yield_FW = yieldFW, W_fruit = Wfruit, Pg, node, harv_n, harv_age,
            W_leaf = Wl, W_stem = Ws, W_root = Wr, yield_DW = Ydw, leaves = Lv, state = s)
end

"""
    load_measured_climate_par(gh_csv, meteo_csv; transmission, lamp_umol, day0_serial, dt) -> MeasuredSeason

Note 69. Canopy light from the organisers' formula (AGC-2018 'Autonomous greenhouse control', p. 8):
    PAR_in = PARout * transmission * (1 - 0.98 BlackScr) * (1 - 0.25 EnScr) + lamp_umol * lamps     [umol m-2 s-1]
converted to W m-2 PAR with 4.6 umol/J (the photosynthesis convention). PARout missing -> 2.16 * Iglob.
transmission: greenhouse light transmission (organisers: about 0.60-0.70); lamp_umol: installed lamp PAR (187).
"""
function load_measured_climate_par(gh_csv::AbstractString, meteo_csv::AbstractString;
                                   transmission::Float64 = 0.65, lamp_umol::Float64 = 187.0,
                                   day0_serial::Float64 = 43326.0, dt::Float64 = 300.0, screen_p = nothing)
    # Note 101 (default nothing = organisers' formula): screen_p = climate parameters with the scr_* optics -> the
    # screen factor is screen_factors(EnScr, BlackScr, screen_p).P.S, the same PAR factor the climate model uses.
    gh = DataFrame(CSV.File(gh_csv)); mt = DataFrame(CSV.File(meteo_csv)); n = min(nrow(gh), nrow(mt))
    isnum(x) = (x isa Number) && !isnan(x)
    v(x, dflt) = isnum(x) ? Float64(x) : dflt
    buckets = Dict{Int, Vector{NTuple{3,Float64}}}()
    for i in 1:n
        g = gh.GHtime[i]; isnum(g) || continue
        Tc = gh.Tair[i]; Co = gh.CO2air[i]; (isnum(Tc) && isnum(Co)) || continue
        d = _excel_day(g, day0_serial); d >= 0 || continue
        Ig  = v(mt.Iglob[i], 0.0); Pout = v(mt.PARout[i], 2.16 * Ig)
        bs  = clamp(v(gh.BlackScr[i], 0.0) / 100, 0.0, 1.0); es = clamp(v(gh.EnScr[i], 0.0) / 100, 0.0, 1.0)
        lam = clamp(v(gh.AssimLight[i], 0.0) / 100, 0.0, 1.0)
        sf  = screen_p === nothing ? (1 - 0.98 * bs) * (1 - 0.25 * es) : screen_factors(es, bs, screen_p).P.S
        par = max(Pout, 0.0) * transmission * sf + lamp_umol * lam
        push!(get!(buckets, d, NTuple{3,Float64}[]), (Tc + 273.15, _co2_ppm_to_mgm3(Co, Tc), par / 4.6))
    end
    ndays = isempty(buckets) ? 0 : maximum(keys(buckets)) + 1
    day_samples = Vector{Matrix{Float64}}(undef, ndays)
    for d in 0:(ndays - 1)
        vv = get(buckets, d, NTuple{3,Float64}[]); M = Matrix{Float64}(undef, 3, length(vv))
        for (j, t) in enumerate(vv); M[1, j] = t[1]; M[2, j] = t[2]; M[3, j] = t[3]; end
        day_samples[d + 1] = M
    end
    return MeasuredSeason(ndays, day_samples, dt)
end

"""
    load_fruit_keep(cm_csv; ndays, day0_date) -> Vector{Float64}

Note 69. Grower policy 'fruits maintained per internode' (CropManagement.csv column `Pruning`, %), as a per-day share
(index d+1), held constant within the week and carried forward over missing weeks (first value before the first week).
"""
function load_fruit_keep(cm_csv::AbstractString; ndays::Int = 116, day0_date::Date = Date(2018, 8, 14),
                         as_removal::Bool = false)
    # Note 97: `Pruning` is the % of fruits REMOVED (Hemming et al. 2019: "removing alternately 50% and 67% of the
    # fruits"); as_removal = true returns the kept share 1 - Pruning/100. Default false reproduces Notes 69-95.
    df = DataFrame(CSV.File(cm_csv))
    wk = [Dates.value((Date(2018, 1, 1) + Day((Int(w) - 1) * 7)) - day0_date) for w in df.weeks]
    pv = [x isa Number && !isnan(x) ? (as_removal ? 1.0 - Float64(x) / 100 : Float64(x) / 100) : NaN for x in df.Pruning]
    first_ok = findfirst(!isnan, pv); first_ok === nothing && error("no Pruning values")
    keep = fill(pv[first_ok], ndays); cur = pv[first_ok]
    for d in 0:(ndays - 1)
        i = findlast(<=(d), wk)
        if i !== nothing && !isnan(pv[i]); cur = pv[i]; end
        keep[d + 1] = cur
    end
    return keep
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
