# =============================================================================
# crop_measured.jl -- the crop driven by measured climate: simulation, observations, likelihood, climate inputs
# =============================================================================
# Port of V2.2 measured.jl (canopy path), test/calibrate_canopy_growth.jl (PARBELOW=1, tag n105) and
# test/crop_posterior_series.jl. Data (AGC-2018, per team): Greenhouse_climate.csv, meteo.csv, Production.csv,
# CropManagement.csv. Day 0 = 14 Aug 2018 (Excel serial 43326).

"Season of measured climate by day; each day a 3 x N matrix [air T (K); CO2 (mg m-3); canopy-top PAR (W m-2)]."
struct MeasuredSeason
    ndays::Int
    day_samples::Vector{Matrix{Float64}}
    dt::Float64                    # seconds between samples
end

_excel_day(serial, day0) = Int(floor(serial - day0))
_co2_ppm_to_mgm3(ppm, Tc) = ppm * 536.4 / (Tc + 273.15)

"""
    load_measured_climate_par(gh_csv, meteo_csv; transmission, lamp_umol) -> MeasuredSeason

Canopy-top PAR from the organisers' formula:
    PAR = PARout x transmission x (1 - 0.98 BlackScr)(1 - 0.25 EnScr) + lamp_umol x lamps   [umol m-2 s-1] / 4.6
(PARout missing -> 2.16 Iglob). Rows need air T and CO2.
"""
function load_measured_climate_par(gh_csv::AbstractString, meteo_csv::AbstractString;
                                   transmission::Float64 = 0.65, lamp_umol::Float64 = 187.0,
                                   day0_serial::Float64 = 43326.0, dt::Float64 = 300.0)
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
        sf  = (1 - 0.98 * bs) * (1 - 0.25 * es)
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
    daily_assimilation_measured(M, lai, p, dt, k_ext) -> Pg [g CH2O m-2 d-1]

FvCB on each sample with the PAR on both canopy faces (measured top PAR x par_both_faces_factor(lai)) and the
effective leaf area (1 - exp(-k_ext lai)) / k_ext; constant leaf dark respiration.
"""
function daily_assimilation_measured(M::Matrix{Float64}, lai::Float64, p::NamedTuple, dt::Float64, k_ext::Float64)
    size(M, 2) == 0 && return 0.0
    s = 0.0
    g = par_both_faces_factor(lai, p)
    @inbounds for j in 1:size(M, 2)
        s += assimilation(M[1, j], M[2, j], g * M[3, j], lai_eff(lai, k_ext), p)
    end
    return s * dt * 1e-3 * (30.0 / 44.0)
end

"Daily mean air temperature [C] (20 C if no samples)."
daily_mean_T(M::Matrix{Float64}) = size(M, 2) == 0 ? 20.0 : sum(@view M[1, :]) / size(M, 2) - 273.15

"""
    load_fruit_keep(cm_csv; ndays) -> Vector{Float64}

Grower's kept-fruit share per day (index d+1) = 1 - Pruning/100 (CropManagement `Pruning` = % of fruits removed),
held within the week and carried over missing weeks.
"""
function load_fruit_keep(cm_csv::AbstractString; ndays::Int = 116, day0_date::Date = Date(2018, 8, 14))
    df = DataFrame(CSV.File(cm_csv))
    wk = [Dates.value((Date(2018, 1, 1) + Day((Int(w) - 1) * 7)) - day0_date) for w in df.weeks]
    pv = [x isa Number && !isnan(x) ? 1.0 - Float64(x) / 100 : NaN for x in df.Pruning]
    first_ok = findfirst(!isnan, pv); first_ok === nothing && error("no Pruning values")
    keep = fill(pv[first_ok], ndays); cur = pv[first_ok]
    for d in 0:(ndays - 1)
        i = findlast(<=(d), wk)
        if i !== nothing && !isnan(pv[i]); cur = pv[i]; end
        keep[d + 1] = cur
    end
    return keep
end

"(days, N_leaves) from CropManagement.csv: cumulative leaves per stem at the week starts."
function load_cropmanagement(csv_path::AbstractString; day0_date::Date = Date(2018, 8, 14))
    df = DataFrame(CSV.File(csv_path))
    days    = [Dates.value((Date(2018, 1, 1) + Day((Int(w) - 1) * 7)) - day0_date) for w in df.weeks]
    Nleaves = Float64.(df.N_leaves)
    keep    = .!isnan.(Nleaves) .& (days .>= 0)
    return days[keep], Nleaves[keep]
end

"(days, cumulative yield kg FW m-2) from Production.csv."
function load_production(prod_csv::AbstractString; day0_serial::Float64 = 43326.0)
    df = DataFrame(CSV.File(prod_csv))
    return _excel_day.(Float64.(df.time), day0_serial), Float64.(df.Total_Prod_cum)
end

"""
    simulate_crop_measured(p, gp, G, season, keep) -> NamedTuple of daily series

p: photosynthesis parameters (climate_n79.json with the crop J_max); gp: crop_n105.json values; G: cucumber_generic.
"""
function simulate_crop_measured(p::NamedTuple, gp, G::CucumberGeneric, season::MeasuredSeason, keep)
    s = init_crop_state(gp)
    nd = season.ndays
    day = Int[]; LAI = Float64[]; yieldFW = Float64[]; Wfruit = Float64[]; Pg = Float64[]; node = Float64[]
    harv_n = Float64[]; harv_age = Float64[]; Wl = Float64[]; Ws = Float64[]; Wr = Float64[]; Ydw = Float64[]
    Lv = Float64[]
    for d in 0:(nd - 1)
        M = season.day_samples[d + 1]
        Pgd = daily_assimilation_measured(M, s.LAI, p, season.dt, gp.k_ext)
        parMJ = size(M, 2) == 0 ? NaN : sum(@view M[3, :]) * season.dt * 1e-6
        kp = Float64(keep[min(d + 1, length(keep))])
        hn, ha = crop_step!(s, Pgd, daily_mean_T(M), parMJ, kp, d, gp, G)
        push!(day, d); push!(LAI, s.LAI); push!(yieldFW, s.yield_FW)
        push!(Wfruit, isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits))
        push!(Pg, Pgd); push!(node, s.node); push!(harv_n, hn); push!(harv_age, hn > 0 ? ha / hn : NaN)
        push!(Wl, s.W_leaf); push!(Ws, s.W_stem); push!(Wr, s.W_root); push!(Ydw, s.yield_DW); push!(Lv, s.leaves)
    end
    return (; day, LAI, yield_FW = yieldFW, W_fruit = Wfruit, Pg, node, harv_n, harv_age,
            W_leaf = Wl, W_stem = Ws, W_root = Wr, yield_DW = Ydw, leaves = Lv, state = s)
end

# --- observations and likelihood (calibrate_canopy_growth.jl) --------------------------------------------------------
const CROP_QOI = ["W_pot", "E_fruit", "app_beta1", "V20", "k_m_leaf", "k_m_stem", "k_m_fruit", "gen_leaf_target",
                  "gen_N0", "set_r_low", "set_r_high", "DMC_fruit", "J_max"]
const CROP_NOISE = ["sigma_y", "sigma_leaf", "sigma_fd", "sigma_n"]

_cum_at(dd, cum, d) = (d <= dd[1] ? float(cum[1]) : d >= dd[end] ? float(cum[end]) :
    let i = searchsortedlast(dd, d); dd[i] == d ? float(cum[i]) : cum[i] + (d - dd[i]) / (dd[i+1] - dd[i]) * (cum[i+1] - cum[i]) end)
_weekly_incr(dd, cum, e) = diff([_cum_at(dd, cum, x) for x in e])
_at_day(v, d) = v[min(d, length(v) - 1) + 1]
function _week_age(r, d0)
    n = 0.0; a = 0.0
    for d in d0:min(d0 + 6, length(r.day) - 1)
        hn = r.harv_n[d + 1]; hn > 0 && (n += hn; a += hn * r.harv_age[d + 1])
    end
    n > 0 ? a / n : NaN
end

"""
    crop_observations(inv; team) -> NamedTuple

Measured season (organisers' PAR formula, transmission and lamp from crop_n105.json defaults), kept-fruit share and
the four observables: weekly yield, weekly harvested fruit number, leaves per stem, fruit development time.
"""
function crop_observations(inv::AbstractString; team::AbstractString = "AiCU", transmission = 0.65, lamp_umol = 187.0)
    data = joinpath(inv, "invernaderos2.0", "data")
    season = load_measured_climate_par(joinpath(data, team, "Greenhouse_climate.csv"), joinpath(data, "meteo.csv");
                                       transmission = transmission, lamp_umol = lamp_umol)
    cm_csv = joinpath(data, team, "CropManagement.csv")
    keep = load_fruit_keep(cm_csv; ndays = season.ndays)
    days_obs, y_obs = load_production(joinpath(data, team, "Production.csv"))
    last = min(season.ndays - 1, days_obs[end])
    edges = let e = collect(0:7:last); e[end] < last && push!(e, last); e end
    dy = _weekly_incr(days_obs, y_obs, edges)
    prod = DataFrame(CSV.File(joinpath(data, team, "Production.csv")))
    numv(x) = x isa Number && !isnan(x) ? Float64(x) : 0.0
    cnt_days = [Int(round(t - 43326)) for t in prod.time]
    cnt_cum  = cumsum([numv(a) + numv(b) for (a, b) in zip(prod.ProdA_num, prod.ProdB_num)])
    dn = _weekly_incr(cnt_days, cnt_cum, edges)
    cm_days, cm_N = load_cropmanagement(cm_csv)
    k = cm_days .<= season.ndays - 1; cm_days = cm_days[k]; cm_N = cm_N[k]
    cmdf = DataFrame(CSV.File(cm_csv))
    fd_start = [Dates.value((Date(2018, 1, 1) + Day((Int(w) - 1) * 7)) - Date(2018, 8, 14)) for w in cmdf.weeks]
    fd_val   = [x isa Number && !isnan(x) ? Float64(x) : NaN for x in cmdf.FruitGrw]
    okfd     = .!isnan.(fd_val) .& (fd_start .>= 0) .& (fd_start .+ 6 .<= season.ndays - 1)
    return (; season, keep, edges, days_obs, y_obs, dy, dn, cm_days, cm_N, fd_start = fd_start[okfd], fd_val = fd_val[okfd],
              last)
end

"Residual vectors (weekly yield, leaves, fruit development time, weekly fruit number) of a simulation r."
function crop_residuals(r, O)
    ry = O.dy .- _weekly_incr(r.day, r.yield_FW, O.edges)
    rl = O.cm_N .- [_at_day(r.node, d) for d in O.cm_days]
    rf = [isnan(_week_age(r, d)) ? v : v - _week_age(r, d) for (d, v) in zip(O.fd_start, O.fd_val)]
    rn = O.dn .- _weekly_incr(r.day, cumsum(r.harv_n), O.edges)
    return ry, rl, rf, rn
end

"Run the crop at z = values of CROP_QOI (pc = photosynthesis parameters, gp0 = crop_n105.json values)."
function crop_run(z, pc, gp0, O)
    gp = update_params(gp0, CROP_QOI[1:end-1], z[1:length(CROP_QOI)-1])
    G  = cucumber_generic(gp.W_pot, gp.E_fruit, gp.app_beta1, gp.V20, gp)
    p  = update_params(pc, ["J_max"], [z[length(CROP_QOI)]])
    return simulate_crop_measured(p, gp, G, O.season, O.keep)
end

"Crop energy at z = values of CROP_QOI followed by sigma_y, sigma_leaf, sigma_fd, sigma_n."
function crop_energy(z, pc, gp0, O)
    r = try crop_run(z, pc, gp0, O) catch; return Inf end
    ry, rl, rf, rn = crop_residuals(r, O); nq = length(CROP_QOI)
    σy, σl, σf, σn = z[nq+1], z[nq+2], z[nq+3], z[nq+4]
    length(ry) * log(σy) + sum(abs2, ry) / (2σy^2) + length(rl) * log(σl) + sum(abs2, rl) / (2σl^2) +
    length(rf) * log(σf) + sum(abs2, rf) / (2σf^2) + length(rn) * log(σn) + sum(abs2, rn) / (2σn^2)
end

"""
    crop_series(r, gp) -> DataFrame

Daily climate inputs from a crop simulation (per m2 floor): LAI_floor, CropResp20 = non-leaf maintenance + growth
respiration at 20 C [mg CO2 m-2 s-1] (leaf dark respiration is in the FvCB model), FW_mass_kg = standing fresh mass
of stems, roots and fruits on the plant [kg m-2].
"""
function crop_series(r, gp)
    DM   = r.W_leaf .+ r.W_stem .+ r.W_root .+ r.W_fruit .+ r.yield_DW
    dDM  = max.(vcat(0.0, diff(DM)), 0.0)
    Rm20 = 2.0^((20 - 25) / 10) .* (gp.k_m_stem .* r.W_stem .+ gp.k_m_root .* r.W_root .+ gp.k_m_fruit .* r.W_fruit)
    Rg   = (gp.asrq - 1) .* dDM
    conv = 44 / 30 * 1000 / 86400
    CR20 = 1.0 .* (Rm20 .+ Rg) .* conv
    return DataFrame(day = r.day, LAI_floor = 1.0 .* r.LAI, CropResp20 = CR20, Rm20_gCH2O = Rm20, Rg_gCH2O = Rg,
                     Pg_gCH2O = r.Pg, yield_FW = r.yield_FW, node = r.node, W_leaf = r.W_leaf, W_stem = r.W_stem,
                     W_root = r.W_root, W_fruit = r.W_fruit, yield_DW = r.yield_DW, harv_n = r.harv_n, LAI = r.LAI,
                     FW_mass_kg = 1.0 .* (r.W_stem ./ gp.dmc_stem .+ r.W_root ./ gp.dmc_root .+ r.W_fruit ./ gp.DMC_fruit) ./ 1000)
end
