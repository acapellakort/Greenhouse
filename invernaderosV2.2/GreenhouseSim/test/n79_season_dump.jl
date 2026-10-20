# =============================================================================
# n79_season_dump.jl -- full-season simulation of the n79 set-up at the PRIOR CENTRE (Note 82; diagnostic only)
# =============================================================================
#     julia --project=. test/crop_posterior_series.jl        (re-run once: adds node / organ weights to the crop CSV)
#     julia --project=. test/n79_season_dump.jl              (~1-2 min)
# Writes test/figures/n79_season_centre.csv (5-min states, observations, drivers, cover T, condensation and
# transpiration fluxes) for two variants: "on" (cover condensation + crop respiration, as n79) and "off" (both off),
# and test/figures/n79_daily_water.csv (daily transpiration, condensation, observed water use on both area bases).
# =============================================================================
ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n79"; haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

const GS = GreenhouseSim
function p1_flux(T1, T2, C1, V1, I9, LAI, p)        # kg m-2 s-1 (climate_rhs, stomata_fix form)
    (; gamma, gamma2, gamma3, gamma4, gamma5, delta1, delta2, delta3, delta4, delta5, delta6, delta7, rho3, alpha5) = p
    q2  = Pws(T1); q7 = 1 / (1 + exp((I9 - delta1) / gamma5)); q10 = (I9 + delta2) / (I9 + delta3)
    q4  = 1 + (delta5 + (delta4 - delta5) * q7) * (C1 * T2 / 536.4 - 200)^2
    q5  = 1 + (delta7 + (delta6 - delta7) * q7) * (q2 - V1)^2
    q1  = 2 * rho3 * alpha5 * LAI / (gamma * gamma2 * (gamma3 + gamma4 * q10 * q4 * q5))
    return q1 * (q2 - V1)
end
Ig = A(Ig_v).(predtimes); Tsk = A(Tsky_v).(predtimes); Tou = A(Tout_v).(predtimes); LAIt = A(LAI_row).(predtimes)
ESt = A(U1v).(predtimes); BSt = A(Blackv).(predtimes); LMt = A(U12v).(predtimes); U10t = A(U10v).(predtimes)
Leet = A(LeeZ).(predtimes); Wndt = A(WindZ).(predtimes); Tpt = A(Tpipe_meas).(predtimes); CRt = A(CropRespv).(predtimes)
Cobs = A(CO2_obs).(predtimes); Wsp = A(Wind_v).(predtimes)
pC = update_params(params_cl, CQ_MODEL, centre[1:NM])
variants = ("on" => pC, "off" => update_params(pC, ["cover_cond", "co2_resp"], [0.0, 0.0]))
out = DataFrame(day = (predtimes .- t0) ./ 86400, Iglob = Ig, Tout = Tou .- 273.15, Tsky = Tsk .- 273.15, wind = Wsp,
                EnScr = ESt, BlackScr = BSt, lamps = LMt, U10 = U10t, VentLee = Leet, VentWind = Wndt,
                Tpipe = Tpt .- 273.15, LAI_floor = LAIt, CropResp20 = CRt,
                Tobs = Tobs_g .- 273.15, RHobs = RHobs_g, CO2obs = Cobs)
daily = DataFrame(day = 0:(N_DAYS - 1))
trans9 = (1 - pC.eta1) * pC.tau1 * pC.eta2
for (tag, p) in variants
    us = simdays(p, days, u0_cl, RELTOL); us === nothing && error("solve failed ($tag)")
    nm = ["T1", "T2", "C1", "V1", "PIint", "T5", "U7state", "T8", "Tdeep", "Tmass"]
    for k in 1:size(us, 1)
        v = us[k, :]; nm[k] in ("T1", "T2", "T5", "T8", "Tdeep", "Tmass") && (v = v .- 273.15)
        out[!, "$(nm[k])_$tag"] = v
    end
    out[!, "RH_$tag"] = rhf.(us[2, :], us[4, :])
    Tc = similar(Ig); cond = similar(Ig); tr = similar(Ig)
    for c in eachindex(Ig)
        T2 = us[2, c]; V1 = us[4, c]; I4 = Tsk[c]; I5 = Tou[c]; U1 = clamp(ESt[c], 0, 1)
        kR = 4.0 * p.epsil3 * p.sigma * (0.5 * (T2 + I4))^3
        tc = (p.solar_split > 0.5) ? GS.cover_T(T2, I4, I5, Ig[c], p) : (kR * (T2 + I4) + p.h_cover_out * I5) / (2kR + p.h_cover_out)
        p.cover_cond > 0.5 && (tc = GS.cover_latent(tc, T2, V1, U1, kR, p))
        Tc[c] = tc - 273.15
        dv = V1 - Pws(tc); fcov = 1 - U1 * (1 - p.vent_scr_perm)
        cond[c] = dv > 0 ? p.cover_area * fcov * 6.4e-9 * 1.7 * abs(T2 - tc)^0.33 * dv * 3.6e6 : 0.0   # g m-2 h-1 (potential if off)
        Ipar = trans9 * Ig[c] * (1 - 0.98 * clamp(BSt[c], 0, 1)) * (1 - 0.25 * U1) + p.lamp_par_W * clamp(LMt[c], 0, 1)
        tr[c] = p1_flux(us[1, c], T2, us[3, c], V1, Ipar, LAIt[c], p) * 3.6e6                          # g m-2 h-1
    end
    out[!, "Tcover_$tag"] = Tc; out[!, "cond_gh_$tag"] = cond; out[!, "transp_gh_$tag"] = tr
    dd = floor.(Int, out.day); dd[1] = 0
    daily[!, "transp_L_$tag"] = [sum(tr[(dd .== d) .& (1:length(dd) .> 1)]) * 300 / 3600 / 1000 for d in 0:(N_DAYS - 1)]
    daily[!, "cond_L_$tag"]   = [sum(cond[(dd .== d) .& (1:length(dd) .> 1)]) * 300 / 3600 / 1000 for d in 0:(N_DAYS - 1)]
end
irr = DataFrame(CSV.File(joinpath(DATA, TEAM, "Irrigation.csv"))); numx(x) = x isa Number ? Float64(x) : NaN
net = fill(NaN, N_DAYS)
for r in eachrow(irr)
    d = Int(floor(numx(r.time) - DAY0_SERIAL)); 0 <= d < N_DAYS || continue
    net[d + 1] = numx(r.water) - (isnan(numx(r.drain)) ? 0.0 : numx(r.drain))
end
harv = [cum_at(days_obs, y_obs, d) - cum_at(days_obs, y_obs, d - 1) for d in 0:(N_DAYS - 1)]; harv[1] = 0.0
for (lab, basis) in (("floor", 96.0), ("grow", 76.8))
    daily[!, "Eobs_$lab"] = net .* (basis / 96) .- harv .* 0.96 .* (76.8 / 96) .- 0.06 * 76.8 / 96
end
mkpath(joinpath(@__DIR__, "figures"))
CSV.write(joinpath(@__DIR__, "figures", "n79_season_centre.csv"), out)
CSV.write(joinpath(@__DIR__, "figures", "n79_daily_water.csv"), daily)
@printf("\nwrote figures/n79_season_centre.csv (%d rows) and figures/n79_daily_water.csv\n", nrow(out))
@printf("season totals (L/m2): transpiration on %.0f / off %.0f; condensation on %.0f (potential, off: %.0f); water use floor %.0f grow %.0f\n",
        sum(daily.transp_L_on), sum(daily.transp_L_off), sum(daily.cond_L_on), sum(daily.cond_L_off),
        sum(filter(!isnan, daily.Eobs_floor)), sum(filter(!isnan, daily.Eobs_grow)))
