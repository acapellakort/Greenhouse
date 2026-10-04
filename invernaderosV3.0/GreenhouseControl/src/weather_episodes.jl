# =============================================================================
# weather_episodes.jl -- weather-year sampling + noisy forecast for RL episodes
# =============================================================================
# Domain randomization: each episode samples one FULLY-COVERED weather year from
# the 32-year Holland record (the file has gaps -- e.g. 2018 ends mid-July -- so
# incomplete windows are filtered out). The agent additionally sees a NOISY
# multi-day forecast of the drivers that move control (outside T, radiation,
# wind); the twin always runs on the TRUE sampled-year weather.

using Random, Statistics

# ---- weather bank: parse the meteo file ONCE, list covered start-years -------
struct WeatherBank
    df::Any                 # prepared meteo frame (GreenhouseSim.prepare_meteo)
    plant_month::Int
    plant_day::Int
    ndays::Int
    years::Vector{Int}
    sky_from_clouds::Bool
end

"Years whose [plant_md, +ndays] window is fully covered (hourly) AND non-degenerate.
The record contains placeholder years filled entirely with zeros (e.g. 2023: Tout=0,
radiation=0 for the whole window) that pass a row-count check but are unusable; we
reject any window whose radiation is all ~zero or whose temperature is ~constant."
function valid_start_years(df, m::Int, d::Int, ndays::Int)
    need = ndays * 24 - 12           # hourly rows expected, small slack
    ut  = df.Unix_Time
    rad = df.External_Radiation; T = df.Temperature
    out = Int[]
    for y in sort(unique(df.Year))
        u0 = datetime2unix(DateTime(y, m, d)); u1 = u0 + ndays * 86400.0
        mask = (ut .>= u0) .& (ut .<= u1)
        count(mask) >= need || continue
        (maximum(rad[mask]) > 50.0 && std(T[mask]) > 0.5) || continue   # reject all-zero/constant windows
        push!(out, y)
    end
    return out
end

function WeatherBank(csv_path::AbstractString; plant_month = 8, plant_day = 14,
                     ndays = 100, sky_from_clouds = true)
    df  = prepare_meteo(csv_path)
    yrs = valid_start_years(df, plant_month, plant_day, ndays)
    isempty(yrs) && error("WeatherBank: no fully-covered $(ndays)-day windows from $plant_month/$plant_day")
    return WeatherBank(df, plant_month, plant_day, ndays, yrs, sky_from_clouds)
end

"Sample (or force) a covered year -> (weather, start_unix, year). Window is padded for forecast lookahead."
function sample_window(bank::WeatherBank; rng = Random.default_rng(), year = nothing, pad_days = 5)
    y     = year === nothing ? rand(rng, bank.years) : year
    start = DateTime(y, bank.plant_month, bank.plant_day)
    stop  = start + Day(bank.ndays + pad_days)
    w = weather_from_df(bank.df, start, stop; sky_from_clouds = bank.sky_from_clouds)
    return w, datetime2unix(start), y
end

# ---- noisy forecast ----------------------------------------------------------
# lead-dependent error std: temp [C], radiation [fraction], wind [m/s].
_sigT(l) = 1.0  + 0.75 * (l - 1)       # 1.0, 1.75, 2.5
_sigR(l) = 0.15 + 0.10 * (l - 1)       # 0.15, 0.25, 0.35
_sigW(l) = 1.0  + 0.5  * (l - 1)       # 1.0, 1.5, 2.0

struct SeasonForecast
    Tout::Vector{Float64}   # true daily-mean outside T [C]  (index = day+1)
    rad::Vector{Float64}    # true daily-mean radiation [W/m2]
    wind::Vector{Float64}   # true daily-mean wind [m/s]
    zT::Matrix{Float64}; zR::Matrix{Float64}; zW::Matrix{Float64}   # [day, lead] episode noise
    H::Int
    skill::Float64          # 0 = perfect forecast ... 1 = nominal ... large = useless
end

"True daily means of the drivers over `ndays` from `start_unix`."
function season_daily_truth(w, start_unix::Float64, ndays::Int)
    T = Float64[]; R = Float64[]; W = Float64[]
    for d in 0:(ndays - 1)
        ts = start_unix .+ d * 86400.0 .+ (0:23) .* 3600.0
        push!(T, mean(w.Tout(t) - 273.15 for t in ts))
        push!(R, mean(w.Iglobal(t)       for t in ts))
        push!(W, mean(w.WindSpeed(t)      for t in ts))
    end
    return T, R, W
end

function SeasonForecast(w, start_unix::Float64, ndays::Int; H = 3, skill = 1.0,
                        rng = Random.default_rng(), pad = 5)
    n = ndays + pad
    T, R, W = season_daily_truth(w, start_unix, n)
    SeasonForecast(T, R, W, randn(rng, n, H), randn(rng, n, H), randn(rng, n, H), H, skill)
end

"Flat noisy forecast row at `day` (0-based): for lead 1..H -> (Tout[C], rad[W/m2], wind[m/s])."
function forecast_row(fc::SeasonForecast, day::Int)
    out = Float64[]
    for l in 1:fc.H
        i  = clamp(day + l + 1, 1, length(fc.Tout))       # target day index
        Tf = fc.Tout[i] + fc.zT[i, l] * _sigT(l) * fc.skill
        Rf = max(fc.rad[i]  * (1 + fc.zR[i, l] * _sigR(l) * fc.skill), 0.0)
        Wf = max(fc.wind[i] +      fc.zW[i, l] * _sigW(l) * fc.skill,  0.0)
        push!(out, Tf, Rf, Wf)
    end
    return out
end
