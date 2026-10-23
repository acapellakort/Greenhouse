# =============================================================================
# export_n106_sim.jl  --  n106 posterior-mean season traces -> JSON
# =============================================================================
# Run from invernaderosN79/GreenhouseTwin/:
#   CROPTAG=n105 CLIMTAG=n106 SEED=1 RELTOL=1e-4 \
#       julia --project=. test/export_n106_sim.jl 2>&1 | tee export_n106.log
# Output: test/posteriors/n106_season_trace_seed<SEED>.json
# =============================================================================
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin, JLD2, Statistics
const GT = GreenhouseTwin

const INV     = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CROPTAG = get(ENV, "CROPTAG", "n105")
const CLIMTAG = get(ENV, "CLIMTAG", "n106")
const SEED    = parse(Int, get(ENV, "SEED", "1"))
const RELTOL  = parse(Float64, get(ENV, "RELTOL", "1e-4"))
const POST    = joinpath(@__DIR__, "posteriors")

cs = let a = joinpath(POST, "AiCU_crop_series_$(CROPTAG).csv"),
         b = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "test", "posteriors",
                      "AiCU_crop_series_$(CROPTAG).csv")
    isfile(a) ? a : b
end

fn = joinpath(POST, "AiCU_climate_$(CLIMTAG)_twin_seed$(SEED).jld2")
isfile(fn) || error("$fn not found — run calibrate_climate.jl full first")

d    = JLD2.load(fn)
post = d["chain"][d["burn_in"]+1:end, 1:end-1]
pm   = vec(mean(post, dims=1))
nm   = length(GT.CLIMATE_QOI)

p0   = load_params(joinpath(@__DIR__, "..", "configfiles", "climate_n79.json"))
S    = load_agc_season(INV; crop_series=cs)
p    = update_params(p0, GT.CLIMATE_QOI, pm[1:nm])

println("Running season_eval (CLIMTAG=$CLIMTAG, seed=$SEED, reltol=$RELTOL) ...")
us, tab = season_eval(p, S; rtol=RELTOL)
println("\nSeason table:")
show(stdout, tab; allrows=true, eltypes=false); println()

# ── reproduce season_eval time-grid ───────────────────────────────────────
const NPD   = 289
const TEVAL_pts = 0.0:300.0:86400.0
days  = collect(0.0:(S.ndays - 1))
pt    = vcat(S.t0, [S.t0 + d * 86400.0 + tj for d in days for tj in TEVAL_pts])

println("\nSeason struct fields: ", fieldnames(typeof(S)))
A(v) = StepFun(S.ts, v)

To  = A(S.Tair_obs).(pt)
Ro  = A(S.RH_obs).(pt)
Co  = A(S.CO2_obs).(pt)
Ig  = A(S.Ig).(pt)

# Try common field names for heater; fall back to zeros
flds = fieldnames(typeof(S))
U1_ts  = :U1  in flds ? A(S.U1).(pt)  :
         :U8  in flds ? A(S.U8).(pt)  : fill(0.0, length(pt))
U12_ts = :U12 in flds ? A(S.U12).(pt) : fill(0.0, length(pt))

T2m = us[2, :]
RHm = rhf.(us[2, :], us[4, :])
C1  = us[3, :]
V1  = us[4, :]

# ── hourly averages ────────────────────────────────────────────────────────
fmt_f(x, n) = round(x, digits=n)
null_or(x, n) = isnothing(x) ? "null" : string(fmt_f(x, n))


# ── hourly averages ────────────────────────────────────────────────────────
fmt_f(x, n) = round(x, digits=n)
null_or(x, n) = isnothing(x) ? "null" : string(fmt_f(x, n))

lines = String[]
for (di, day) in enumerate(days)
    base = (di - 1) * NPD + 2
    for h in 0:23
        hcols = [base + j for j in 0:NPD-1
                 if min(floor(Int, j * 300 / 3600), 23) == h && base + j <= length(pt)]
        isempty(hcols) && continue

        okT  = [c for c in hcols if !isnan(To[c])]
        okR  = [c for c in hcols if !isnan(Ro[c])]
        okTR = [c for c in hcols if !isnan(To[c]) && !isnan(Ro[c])]
        okC  = [c for c in hcols if !isnan(Co[c])]
        half = length(hcols) ÷ 2

        T2o_h = length(okT)  >= half ? mean(To[okT])  - 273.15 : nothing
        RHo_h = length(okR)  >= half ? mean(Ro[okR])            : nothing
        VPo_h = length(okTR) >= half ?
                mean(Ro[okTR] ./ 100 .* Pws.(To[okTR])) / 100 : nothing
        C2o_h = length(okC)  >= half ? mean(Co[okC])            : nothing

        push!(lines, """  {"day":$(fmt_f(day+h/24,4)),"T2m":$(fmt_f(mean(T2m[hcols])-273.15,2)),"T2o":$(null_or(T2o_h,2)),"RHm":$(fmt_f(mean(RHm[hcols]),1)),"RHo":$(null_or(RHo_h,1)),"VPm":$(fmt_f(mean(V1[hcols])/100,2)),"VPo":$(null_or(VPo_h,2)),"C2m":$(fmt_f(mean(C1[hcols]),1)),"C2o":$(null_or(C2o_h,1)),"Ig":$(fmt_f(mean(Ig[hcols]),1)),"U1":$(fmt_f(mean(U1_ts[hcols]),3)),"U12":$(fmt_f(mean(U12_ts[hcols]),3))}""")
    end
end

out = joinpath(POST, "n106_season_trace_seed$(SEED).json")
open(out, "w") do io
    println(io, "[")
    for (i, l) in enumerate(lines)
        i < length(lines) ? println(io, l, ",") : println(io, l)
    end
    println(io, "]")
end
println("\nWrote $(length(lines)) records → $out")
