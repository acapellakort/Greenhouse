# =============================================================================
# diagnose_humidity_co2.jl -- humidity, VP and CO2 scorecard of a calibrated posterior (Note 64; forward only)
# =============================================================================
# Run from the GreenhouseSim folder (~1 min per model, no MCMC):
#     PRIORS=n54 julia --project=. test/diagnose_humidity_co2.jl      (closed loop)
#     PRIORS=n61 julia --project=. test/diagnose_humidity_co2.jl      (open loop, logged vents)
#
# Pooled posterior mean (final .jld2, else checkpoint; burn-in 30 000), full season from day 0.
#   1. RH and VP: bias and hourly RMSE by season period, day / night, and by 2-h hour of day (model AIR node,
#      no measurement term).
#   2. Humid hours: hourly bins with measured RH > 90 % -- hits, misses, false alarms, bias.
#   3. CO2 (never calibrated): bias / RMSE in ppm by period, day / night. Nights with lamps OFF have no dosing
#      and no photosynthesis: there CO2 follows respiration and AIR EXCHANGE only, a check of the vent and
#      leakage model that does not involve the temperature sensor.
#   4. Measured daily CO2 dosage (ResourceCalculations.csv) next to the model's dosing hours (U10 > 0).
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n61")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
using DelimitedFiles

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2"); fc = replace(f, ".jld2" => "_ckpt.jld2")
    src = isfile(f) ? f : isfile(fc) ? fc : ""; isempty(src) && continue
    ch = JLD2.load(src)["chain"]; size(ch, 1) > 31000 || continue
    push!(parts, ch[30001:end, 1:NM]); println("seed $s: ", basename(src))
end
isempty(parts) && error("no $(PRIORS) chains in $POSTDIR")
pm = vec(mean(reduce(vcat, parts), dims = 1)); pP = update_params(params_cl, CQ_MODEL, pm)
us = simdays(pP, days, u0_cl, RELTOL); us === nothing && error("solve failed")
T = us[2, :]; RH = rhf.(us[2, :], us[4, :])
ppm(C, Tk) = C .* Tk ./ 536.4
Cm = ppm(us[3, :], T); Co = ppm(A(CO2_obs).(predtimes), Tobs_g)
VP(Tk, rh) = rh ./ 100 .* Pws.(Tk) ./ 100; VPm = VP(T, RH); VPo = VP(Tobs_g, RHobs_g)

okT = .!isnan.(Tobs_g); okT[1] = false; okR = okT .& .!isnan.(RHobs_g); okC = .!isnan.(Co); okC[1] = false
inset(ds) = in.(cd_, Ref(Set(ds)))
PER = [(0, 8), (9, 48), (49, 100), (101, 109), (110, 114)]
hrc = [floor(Int, mod(((c - 2) % NPD) * 300 / 3600, 24.0)) for c in 1:NCOL]; hrc[1] = 0
lon = A(U12v).(predtimes) .> 0.5
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
function hourly(v, o, m)          # hourly-bin means (model, obs) over masked columns with >= 6 samples
    a = Float64[]; b = Float64[]
    for cols in values(HB); k = [c for c in cols if m[c]]; length(k) >= 6 && (push!(a, mean(v[k])); push!(b, mean(o[k]))); end
    a, b
end
bias(v, o, m) = (x = hourly(v, o, m); isempty(x[1]) ? NaN : mean(x[1] .- x[2]))
rmse(v, o, m) = (x = hourly(v, o, m); isempty(x[1]) ? NaN : sqrt(mean(abs2, x[1] .- x[2])))

println("\n", "="^118)
@printf("  %s  posterior mean, full season, %s\n", PRIORS, OLP ? "OPEN loop (logged vents)" : "CLOSED loop (model vents)")
println("  1. RH (%) and VP (hPa): bias / hourly RMSE  (model air node - measured)")
@printf("  %-10s | %15s %15s | %15s %15s | %15s %15s\n", "period", "RH all", "", "RH day", "RH night", "VP day", "VP night")
for (a, b) in vcat(PER, [(0, 114)])
    md = inset(collect(a:b)) .& okR
    @printf("  d%-3d-%-4d | %+6.2f / %5.2f %15s | %+6.2f / %5.2f | %+6.2f / %5.2f | %+6.2f / %5.2f | %+6.2f / %5.2f\n", a, b,
            bias(RH, RHobs_g, md), rmse(RH, RHobs_g, md), "",
            bias(RH, RHobs_g, md .& dayt), rmse(RH, RHobs_g, md .& dayt), bias(RH, RHobs_g, md .& nite), rmse(RH, RHobs_g, md .& nite),
            bias(VPm, VPo, md .& dayt), rmse(VPm, VPo, md .& dayt), bias(VPm, VPo, md .& nite), rmse(VPm, VPo, md .& nite))
end
com = inset(COMMON_TEST) .& okR
@printf("  common test: RH bias %+.2f, RMSE %.2f %%;  VP bias %+.2f, RMSE %.2f hPa\n", bias(RH, RHobs_g, com), rmse(RH, RHobs_g, com),
        bias(VPm, VPo, com), rmse(VPm, VPo, com))
println("\n  RH bias by hour of day (%):   ", join([@sprintf("%02d:%+5.1f", h, bias(RH, RHobs_g, okR .& (hrc .>= h) .& (hrc .< h + 2))) for h in 0:2:22], "  "))
println("  VP bias by hour of day (hPa): ", join([@sprintf("%02d:%+5.2f", h, bias(VPm, VPo, okR .& (hrc .>= h) .& (hrc .< h + 2))) for h in 0:2:22], "  "))

mo, oo = hourly(RH, RHobs_g, okR)
hit = count((oo .> 90) .& (mo .> 90)); miss = count((oo .> 90) .& (mo .<= 90)); fa = count((oo .<= 90) .& (mo .> 90))
@printf("\n  2. Humid hours (hourly RH > 90 %%): measured %d of %d hours; model hits %d, misses %d, false alarms %d; bias in measured-humid hours %+.1f %%\n",
        count(oo .> 90), length(oo), hit, miss, fa, isempty(oo[oo .> 90]) ? NaN : mean(mo[oo .> 90] .- oo[oo .> 90]))

println("\n  3. CO2 (ppm; never calibrated): bias / hourly RMSE      measured mean")
@printf("  %-10s | %17s | %17s | %17s | %9s\n", "period", "day", "night lamps ON", "night lamps OFF", "obs mean")
for (a, b) in vcat(PER, [(0, 114)])
    md = inset(collect(a:b)) .& okC
    @printf("  d%-3d-%-4d | %+7.0f / %6.0f | %+7.0f / %6.0f | %+7.0f / %6.0f | %9.0f\n", a, b,
            bias(Cm, Co, md .& dayt), rmse(Cm, Co, md .& dayt), bias(Cm, Co, md .& nite .& lon), rmse(Cm, Co, md .& nite .& lon),
            bias(Cm, Co, md .& nite .& .!lon), rmse(Cm, Co, md .& nite .& .!lon), mean(Co[md]))
end
println("  CO2 bias by hour of day (ppm): ", join([@sprintf("%02d:%+4.0f", h, bias(Cm, Co, okC .& (hrc .>= h) .& (hrc .< h + 2))) for h in 0:2:22], " "))
# night lamps-off decay test: correlation of hourly (C - C_out) model vs measured
mo, oo = hourly(Cm, Co, okC .& nite .& .!lon)
@printf("  night lamps-OFF hourly CO2: corr(model, obs) %.2f, mean model %.0f vs obs %.0f ppm (%d hours)\n", cor(mo, oo), mean(mo), mean(oo), length(mo))

rc = joinpath(DATA, TEAM, "ResourceCalculations.csv")
if isfile(rc)
    raw, hdr = readdlm(rc, ',', header = true); hdr = vec(String.(hdr)); j = findfirst(==("CO2_dosage"), hdr)
    dose = [x isa Number ? Float64(x) : NaN for x in raw[:, j]]
    u10on = [mean(A(U10v).(predtimes)[inset([d])] .> 0) * 24 for d in 0:114]
    @printf("\n  4. Measured CO2 dosage: mean %.4f (units of ResourceCalculations, per day; %d days).  Model dosing (U10 > 0): %.1f h/day on average\n",
            mean(filter(!isnan, dose)), count(!isnan, dose), mean(u10on))
end
println("""
  Read: RH/VP errors by hour show whether transpiration or ventilation is off; humid-hour misses matter for disease
  risk and humidity venting. Night lamps-OFF CO2 tests the air exchange independently of the T sensor: a model that
  is too high there ventilates too little, too low too much (or respiration is off). Daytime CO2 mixes dosing
  (crude U10 rule) and photosynthesis and is not yet a fair test.""")
