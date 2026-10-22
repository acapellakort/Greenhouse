# =============================================================================
# agc_data.jl -- AGC-2018 inputs and observations for one compartment (default AiCU), n79 conventions
# =============================================================================
# Port of the input construction of V2.2 test/run_joint_inference_optionB.jl (PIPE2 = 1) and
# test/run_climate_calibration_cl.jl (PRIORS = n79). Data: invernaderos2.0/data/<team>/{Greenhouse_climate.csv,
# vip.csv, ResourceCalculations.csv}, invernaderos2.0/data/meteo.csv, and the crop series written by
# test/crop_series.jl (LAI_floor, CropResp20, FW_mass_kg per day).

const EXCEL_UNIX0    = 25569.0
const SIGMA_SB       = 5.670374419e-8       # used only for the sky temperature from the pyrgeometer
const DAY0_SERIAL    = 43326.0              # 14 Aug 2018 (Excel serial), season day 0
const PIPE_MIN_VALID = 5.0                  # C: logged pipe below this = no reading
const PIPE_DELTA_K   = 25.0
const PIPE_WARM_TOUT = 22.0
const PIPE2_R        = 0.6                  # crop (grow) pipe heat per K^1.32 relative to the rail pipe
const CO2_AMB_PPM    = 390.0
const CO2_THR_PPM    = 50.0
const CO2_INJECT_FRAC = 0.3
const CO2_VALVE_KG_M2_H = 0.0128            # valve capacity, kg m-2 h-1 (dosage -> valve hours)
const JAN1_2018_UNIX = 1_514_764_800.0      # 2018-01-01 00:00 UTC, seconds

"""
    co2_bleiswijk_ppm(t_unix; A_seas=9.0, A_diur=10.0) -> ppm

Outdoor CO₂ model for Bleiswijk 2018 (AGC-2018 site, 52.0°N 4.5°E).
Base 408 ppm (Mauna Loa 2018 annual mean).
Seasonal: ±A_seas ppm, peaks at DOY 120 (≈ May 1, NH maximum).
Diurnal:  ±A_diur ppm, peaks at 05:00 local (CET = UTC+1).
  Amplitude set conservatively at ±10 ppm (outdoor weather station level;
  literature range 5–18 ppm at surface near vegetation).
"""
function co2_bleiswijk_ppm(t_unix::Float64; A_seas::Float64 = 9.0, A_diur::Float64 = 10.0)
    doy  = (t_unix - JAN1_2018_UNIX) / 86400.0 + 1.0   # continuous DOY (1.0 = Jan 1 00:00 UTC)
    hr_local = mod((t_unix + 3600.0) / 3600.0, 24.0)   # UTC+1 (CET); good enough over Aug-Dec window
    seasonal = A_seas * cos(2π * (doy  - 120.0) / 365.0)
    diurnal  = A_diur * cos(2π * (hr_local - 5.0) / 24.0)
    return 408.0 + seasonal + diurnal
end

_serial2unix(s) = (s - EXCEL_UNIX0) * 86400.0
_num(x) = (x isa Number && !(x isa Missing) && !isnan(Float64(x))) ? Float64(x) : NaN
_q2(T) = T > 0 ? 611.21 * exp((18.678 - T / 234.5) * (T / (257.14 + T))) :
                 611.21 * exp((23.036 - T / 333.7) * (T / (279.82 + T)))
_vpout(Tc, RH) = _q2(Tc) * clamp(RH, 1.0, 100.0) / 100.0

"""
    load_agc_season(inv_root; team = "AiCU", ndays = 115, crop_series = <path>) -> NamedTuple

Builds the model inputs (`weather`, `controls`, `t0`) and the observations of one compartment for days
0..ndays-1 after 14 Aug 2018. `crop_series` is the CSV written by test/crop_series.jl.
Returned arrays (one value per logged row, 5 min): ts (Unix), Ig, Tsky, Tout, Wind, VPout, CO2out, U1 (energy screen),
Black, Lee, Wind (vent openings), U10 (CO2 valve), U11 (heating setpoint, K), U12 (lamps), Tpipe (K), LAI, CropResp,
CropFW, and the observations Tair_obs [K], RH_obs [%], CO2_obs [mg m-3]; u0 = initial state for day 0.
CO2out uses the Bleiswijk 2018 seasonal+diurnal model (co2_bleiswijk_ppm); no measured outside CO2 in the dataset.
"""
function load_agc_season(inv::AbstractString; team::AbstractString = "AiCU", ndays::Int = 115,
                         crop_series::AbstractString = "")
    data = joinpath(inv, "invernaderos2.0", "data")
    gh  = DataFrame(CSV.File(joinpath(data, team, "Greenhouse_climate.csv")))
    mt  = DataFrame(CSV.File(joinpath(data, "meteo.csv")))
    vip = DataFrame(CSV.File(joinpath(data, team, "vip.csv")))
    nmv = names(vip)
    vcols(base, zs) = (z = [Symbol("$(base)_$(k)") for k in zs if "$(base)_$(k)" in nmv];
                       isempty(z) ? (base in nmv ? [Symbol(base)] : Symbol[]) : z)
    vc_heat = vcols("HeatTemp_Vip", (1, 2, 3, 5, 6))     # AiCU's own zone (_4) has no heating setpoint
    n = min(nrow(gh), nrow(mt), nrow(vip))
    win_lo = DAY0_SERIAL; win_hi = DAY0_SERIAL + ndays

    ts = Float64[]; Ig = Float64[]; Tsky = Float64[]; Tout = Float64[]; Wind = Float64[]; VPout = Float64[]
    U1 = Float64[]; U10 = Float64[]; U11 = Float64[]; U12 = Float64[]; Tpipe = Float64[]; Black = Float64[]
    Lee = Float64[]; Wnd = Float64[]; gday = Int[]
    Tair_obs = Float64[]; CO2_obs = Float64[]; RH_obs = Float64[]
    for i in 1:n
        g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
        To = _num(mt.Tout[i]); Igr = _num(mt.Iglob[i]); Ws = _num(mt.Windsp[i])
        Rho = _num(mt.Rhout[i]); Py = _num(mt.Pyrgeo[i])
        ToK = (isnan(To) ? 15.0 : To) + 273.15
        Igv = isnan(Igr) ? 0.0 : max(Igr, 0.0)
        Ldn = SIGMA_SB * ToK^4 + (isnan(Py) ? 0.0 : Py)
        Tskv = isnan(Py) ? (-0.4 + 273.15) : clamp((Ldn / SIGMA_SB)^0.25, ToK - 30.0, ToK + 2.0)
        vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i]); en = _num(gh.EnScr[i])
        al = _num(gh.AssimLight[i]); pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
        hv = filter(!isnan, [_num(vip[i, c]) for c in vc_heat])
        Th = isempty(hv) ? 20.0 : mean(hv)
        # heating pipe: logged circuits (both warm -> one equivalent pipe, h ~ dT^1.32); no reading -> Tout + 2 K
        pg_c = (!isnan(pg) && pg > PIPE_MIN_VALID) ? pg : -Inf
        pl_c = (!isnan(pl) && pl > PIPE_MIN_VALID) ? pl : -Inf
        pmax = max(pg_c, pl_c)
        if pg_c > -Inf && pl_c > -Inf
            Tal = _num(gh.Tair[i])
            if !isnan(Tal) && pg_c > Tal && pl_c > Tal
                pmax = Tal + ((pl_c - Tal)^1.32 + PIPE2_R * (pg_c - Tal)^1.32)^(1 / 1.32)
            end
        end
        To_val = isnan(To) ? 15.0 : To
        Tp = pmax > -Inf ? pmax + 273.15 :
             (To_val < Th - 2.0 && To_val < PIPE_WARM_TOUT) ? ToK + 2.0 :     # V2.2: assumed-warm rows -> Tout + 2 K
             To_val + 2.0 + 273.15
        co2vip = _num(vip.CO2_Vip[i])
        Ta = _num(gh.Tair[i]); Co = _num(gh.CO2air[i]); Rh = _num(gh.RHair[i]); bs = _num(gh.BlackScr[i])
        push!(ts, _serial2unix(g)); push!(gday, floor(Int, g))
        push!(Ig, Igv); push!(Tout, ToK); push!(Wind, isnan(Ws) ? 2.0 : Ws)
        push!(VPout, _vpout(To_val, isnan(Rho) ? 80.0 : Rho)); push!(Tsky, Tskv)
        push!(U1, isnan(en) ? 0.0 : clamp(en / 100.0, 0.0, 1.0))
        push!(U12, isnan(al) ? 0.0 : clamp(al / 100.0, 0.0, 1.0))
        push!(U10, (isnan(co2vip) || co2vip <= CO2_AMB_PPM + CO2_THR_PPM) ? 0.0 : CO2_INJECT_FRAC)
        push!(U11, Th + 273.15); push!(Tpipe, Tp)
        push!(Black, isnan(bs) ? 0.0 : clamp(bs / 100, 0.0, 1.0))
        push!(Lee, isnan(vl) ? 0.0 : clamp(vl / 100, 0.0, 1.0)); push!(Wnd, isnan(vw) ? 0.0 : clamp(vw / 100, 0.0, 1.0))
        push!(Tair_obs, isnan(Ta) ? NaN : Ta + 273.15)
        push!(CO2_obs, isnan(Co) ? NaN : Co * 536.4 / ((isnan(Ta) ? 20.0 : Ta) + 273.15))
        push!(RH_obs, isnan(Rh) ? NaN : Rh)
    end
    median(Wind) > 25.0 && (Wind ./= 3.6)               # km/h -> m/s if needed
    t0 = ts[1]
    abs(t0 - _serial2unix(DAY0_SERIAL)) < 1.0 || error("first window row is not at day-0 midnight")
    trel = ts .- t0

    # CO2 valve: on-times from the setpoint rule, scaled per day so that valve hours = logged dosage / capacity
    rc  = DataFrame(CSV.File(joinpath(data, team, "ResourceCalculations.csv")))
    dos = Dict(floor(Int, rc.time[i]) => rc.CO2_dosage[i] for i in 1:nrow(rc)
               if !ismissing(rc.CO2_dosage[i]) && !isnan(rc.CO2_dosage[i]))
    dtr = length(ts) > 1 ? median(diff(ts)) / 3600 : 300 / 3600
    for gd in unique(gday)
        rows = findall(==(gd), gday); hlog = get(dos, gd, NaN) / CO2_VALVE_KG_M2_H
        on = [i for i in rows if U10[i] > 0]; hon = length(on) * dtr
        if isnan(hlog)
            for i in on; U10[i] = 1.0; end
        elseif hon > 0
            sc = hlog / hon; for i in on; U10[i] = min(sc, 1.0); end
        else
            dayr = [i for i in rows if Ig[i] > 50.0]
            for i in dayr; U10[i] = isempty(dayr) ? 0.0 : min(hlog / (length(dayr) * dtr), 1.0); end
        end
    end

    # crop inputs (offline coupling): daily LAI, respiration at 20 C, standing fresh mass
    isfile(crop_series) || error("missing crop series $crop_series -- run test/crop_series.jl first")
    cs = DataFrame(CSV.File(crop_series))
    cday(x) = clamp(floor(Int, x), cs.day[1], cs.day[end]) - cs.day[1] + 1
    LAI      = [cs.LAI_floor[cday((t - t0) / 86400.0)] for t in ts]
    CropResp = [cs.CropResp20[cday((t - t0) / 86400.0)] for t in ts]
    CropFW   = [cs.FW_mass_kg[cday((t - t0) / 86400.0)] for t in ts]
    _conv    = 44.0/30.0 * 1000.0/86400.0                              # g CH2O m-2 d-1 → mg CO2 m-2 s-1
    CropRm20 = [cs.Rm20_gCH2O[cday((t - t0) / 86400.0)] * _conv for t in ts]  # maintenance at 20 C
    CropRg   = [cs.Rg_gCH2O[cday((t - t0) / 86400.0)]   * _conv for t in ts]  # growth respiration

    # Outside CO2: Bleiswijk 2018 seasonal+diurnal model, converted to mg m-3 at measured Tout
    CO2out = [co2_bleiswijk_ppm(ts[j]) * 536.4 / Tout[j] for j in eachindex(ts)]

    A(v) = StepFun(ts, v); Z(v) = StepFun(trel, v); z0 = Z(zeros(length(ts)))
    weather  = WeatherInputs(A(Ig), A(Tsky), A(Tout), A(Wind), A(VPout), A(CO2out), ConstFun(273.15), A(LAI))
    controls = ControlInputs(Z(U1), Z(Black), z0, z0, z0, z0, z0, z0, Z(Lee), Z(Wnd), z0, Z(U10), Z(U11), Z(U12),
                             Z(Tpipe), Z(CropResp), Z(CropFW), Z(CropRm20), Z(CropRg))

    # initial state for day 0: measured air T, CO2, RH; floor 18 C; deep soil and mass at the day-0 mean air T
    i0  = findfirst(!isnan, Tair_obs)
    Ta0 = Tair_obs[i0]; Co0 = isnan(CO2_obs[i0]) ? 575.0 : CO2_obs[i0]; Rh0 = isnan(RH_obs[i0]) ? 80.0 : RH_obs[i0]
    Tmean0 = mean(filter(!isnan, A(Tair_obs).(t0 .+ (0.0:300.0:86400.0))))
    u0 = [Ta0, Ta0, Co0, (Rh0 / 100.0) * Pws(Ta0), 0.0, 18.0 + 273.15, Tmean0, Tmean0]
    return (; team, ndays, t0, ts, trel, weather, controls, u0,
              Ig, Tsky, Tout, Wind, VPout, U1, Black, Lee, WindVent = Wnd, U10, U11, U12, Tpipe, LAI, CropResp, CropFW,
              Tair_obs, RH_obs, CO2_obs)
end
