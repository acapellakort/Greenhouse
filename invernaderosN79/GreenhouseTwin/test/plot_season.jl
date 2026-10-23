# =============================================================================
# test/plot_season.jl  --  export full-season traces for eyeballing
# =============================================================================
# Runs the n79 climate model at the prior centre for the whole season and saves
# a CSV with T2, RH, CO2 (model and observed) on the 5-min grid.
#
# Run from invernaderosN79/GreenhouseTwin (the directory that contains Project.toml):
#   CROPTAG=n105 julia --project=. test/plot_season.jl
#
# Output:  test/season_traces.csv
# =============================================================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin, Printf, Statistics

const GT     = GreenhouseTwin
const INV    = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CJ     = joinpath(@__DIR__, "..", "configfiles", "climate_n79.json")
const CROPTAG = get(ENV, "CROPTAG", "n105")
const RELTOL  = parse(Float64, get(ENV, "RELTOL", "1e-3"))

cs_new = joinpath(@__DIR__, "posteriors", "AiCU_crop_series_$(CROPTAG).csv")
cs_old = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "test", "posteriors",
                  "AiCU_crop_series_$(CROPTAG).csv")
const CROPCSV = isfile(cs_new) ? cs_new : cs_old

const p0     = load_params(CJ)
const PR     = priors_of(CJ)
const Q      = vcat(GT.CLIMATE_QOI, GT.CLIMATE_NOISE)
const S      = load_agc_season(INV; crop_series = CROPCSV)
const centre = [PR[q][2] for q in Q]
const p      = update_params(p0, GT.CLIMATE_QOI, centre[1:length(GT.CLIMATE_QOI)])

println("Running full season at the prior centre (reltol=$RELTOL)…")
us, tab = season_eval(p, S; rtol = RELTOL)
println("season_eval complete. Summary:")
show(stdout, tab; allrows = true, eltypes = false); println()

# ------------------------------------------------------------------
# Reconstruct the 5-min time axis and observed arrays  (same logic
# as season_eval in climate_calibration.jl)
# ------------------------------------------------------------------
using GreenhouseTwin: StepFun, rhf, Pws, TEVAL, NPD

A(v) = StepFun(S.ts, v)
pt   = vcat(S.t0, [S.t0 + d * 86400.0 + tj
                   for d in 0.0:(S.ndays - 1) for tj in TEVAL])

To = A(S.Tair_obs).(pt)
Ro = A(S.RH_obs).(pt)
Co = A(S.CO2_obs).(pt)
Ig = A(S.Ig).(pt)
U1 = A(S.U1).(pt)      # energy screen
U12 = A(S.U12).(pt)    # lamps

T2m  = us[2, :]
RHm  = rhf.(us[2, :], us[4, :])
C1m  = us[3, :]
VPm  = RHm ./ 100 .* Pws.(T2m) ./ 100   # hPa
VPo  = Ro  ./ 100 .* Pws.(To)  ./ 100   # hPa

# day index: -1 for the initial point, 0..ndays-1 for the rest
day_idx = vcat(-1, [(c - 2) ÷ NPD for c in 2:length(pt)])
t_day   = (pt .- S.t0) ./ 86400.0       # fractional days since season start

# ------------------------------------------------------------------
# Write CSV  (one row per 5-min point; skip the initial u0 row d=-1)
# ------------------------------------------------------------------
outpath = joinpath(@__DIR__, "season_traces.csv")
open(outpath, "w") do io
    println(io, "t_day,day,T2_mod,T2_obs,RH_mod,RH_obs,VP_mod_hPa,VP_obs_hPa,CO2_mod,CO2_obs,Iglob,U1,U12")
    for i in 2:length(pt)        # skip i=1 (the S.t0 seed point)
        To_i  = isnan(To[i])  ? "" : @sprintf("%.4f", To[i])
        Ro_i  = isnan(Ro[i])  ? "" : @sprintf("%.2f", Ro[i])
        VPo_i = isnan(Ro[i])  ? "" : @sprintf("%.4f", VPo[i])
        Co_i  = isnan(Co[i])  ? "" : @sprintf("%.2f", Co[i])
        @printf(io, "%.5f,%d,%.4f,%s,%.2f,%s,%.4f,%s,%.2f,%s,%.2f,%.3f,%.3f\n",
                t_day[i], day_idx[i],
                T2m[i], To_i,
                RHm[i], Ro_i,
                VPm[i], VPo_i,
                C1m[i], Co_i,
                Ig[i], U1[i], U12[i])
    end
end
println("Wrote $outpath  ($(length(pt)-1) rows)")
