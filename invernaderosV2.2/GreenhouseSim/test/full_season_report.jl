# =============================================================================
# full_season_report.jl -- how good is the climate model over the WHOLE AiCU season? (overview, Notes 33-40)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     julia --project=. test/full_season_report.jl
#
# Simulates days 0-114 once per configuration, each at its pooled posterior mean (seeds 1+2):
#   OL33  Note 33  open loop: logged vents replayed, all Note 32 physics fixes on
#   CL36  Note 36  closed-loop vents (first closed-loop calibration)
#   CL38  Note 38  closed-loop vents, Note 34/37 controller priors   <- current best model
#   NV39  Note 39  CL38 + logged vents at night (what-if: night ventilation right; NOT a forecast model)
# Writes
#   test/figures/AiCU_season_T_daily.png      daily max / mean / min air T, measured vs models
#   test/figures/AiCU_season_bias.png         daily night and day bias, TRAIN/VAL blocks shaded
#   test/figures/AiCU_season_vp_vent.png      daily vapour pressure and vent opening
#   test/figures/AiCU_season_days.png         5-min traces on 8 representative days
#   test/figures/AiCU_season_daily.csv        one row per day (for the overview page)
#   test/figures/AiCU_season_summary.csv      RMSE / bias by set and configuration
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n37"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
using DelimitedFiles

# ── posterior means ──────────────────────────────────────────────────────────
function post_mean(tag)
    ms = Vector{Vector{Float64}}(); names = String[]
    for s in (1, 2)
        f = joinpath(POSTDIR, "$(tag)_seed$(s).jld2"); isfile(f) || error("missing $f")
        d = JLD2.load(f); q = String.(d["QoI"]); names = [x for x in q if !startswith(x, "sigma")]
        push!(ms, vec(mean(d["chain"][d["burn_in"]+1:end, 1:length(names)], dims = 1)))
    end
    names, (ms[1] .+ ms[2]) ./ 2
end
n33, x33 = post_mean("AiCU_climcal_blocks")
n36, x36 = post_mean("AiCU_climcal_cl_blocks")
n38, x38 = post_mean("AiCU_climcal_cl_blocks_n37")

FIX33N = ["nu4", "U8_min_leak", "heat_ctrl", "hc_Kp", "hc_Ki", "vent_screen", "stomata_fix", "vapour_fix", "cond_heat"]
FIX33V = [4.46e-5, 0.0, 1.0, 10.0, 10.0 / 1800, 1.0, 1.0, 1.0, 1.0]
p33 = update_params(update_params(params, FIX33N, FIX33V), n33, x33)
p36 = update_params(params_cl, n36, x36)
p38 = update_params(params_cl, n38, x38)
p39 = update_params(p38, ["vc_night_log"], [1.0])
ctrl_plain = ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas))

function simrun(p, ctrl, u0x)
    us = Matrix{Float64}(undef, length(u0x), 1 + length(days) * NPD); us[:, 1] = u0x
    u = copy(u0x); col = 2
    for d in days
        ctx = GreenhouseSim.SimContext(p, weather_c, ctrl, t0 + 86400.0 * d, 86400.0 * d)
        sol = GreenhouseSim.solve(GreenhouseSim.ODEProblem(GreenhouseSim.rhs!, u, tspan, ctx), GreenhouseSim.Tsit5();
                                  saveat = teval, reltol = RELTOL)
        length(sol.t) == NPD || error("solver stopped on day $d")
        S = Array(sol); us[:, col:col+NPD-1] = S; col += NPD; u = S[:, end]
    end
    us
end
cfgs = [("OL33", "Note 33 open loop",   p33, ctrl_plain, copy(u0), false),
        ("CL36", "Note 36 closed loop", p36, ctrl_x,     u0_8,     true),
        ("CL38", "Note 38 closed loop", p38, ctrl_x,     u0_8,     true),
        ("NV39", "CL38 + logged night vents (what-if)", p39, ctrl_x, u0_8, true)]
R = Dict{String,NamedTuple}()
println("\nRunning $(length(cfgs)) full-season simulations (reltol $(RELTOL)) ...")
for (k, lab, p, ct, u0x, closed) in cfgs
    tt = @elapsed (us_ = simrun(p, ct, u0x))
    R[k] = (T = us_[2, :], RH = rhf.(us_[2, :], us_[4, :]), U = closed ? clamp.(us_[7, :], 0.0, 1.0) : A(U8z).(predtimes))
    @printf("  %-5s %-38s %.1f s\n", k, lab, tt)
end
KEYS = [c[1] for c in cfgs]; LABS = Dict(c[1] => c[2] for c in cfgs)

# ── masks and helpers ────────────────────────────────────────────────────────
okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g); okU = .!isnan.(U8obs); okU[1] = false
VP(T, RH) = RH ./ 100 .* Pws.(T) ./ 100
VPo = VP(Tobs_g, RHobs_g)
setname(d) = d in TRAIN_DAYS ? "TRAIN" : d in VAL_DAYS ? "VAL" : d in COMMON_TEST ? "COMMON" : "TEST"
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
function hourly_rmse(v, o, mask)
    e = Float64[]
    for cols in values(HB); k = [c for c in cols if mask[c]]; length(k) >= 6 && push!(e, mean(v[k]) - mean(o[k])); end
    isempty(e) ? NaN : sqrt(mean(abs2, e))
end
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
rmse_(v, o, m) = any(m) ? sqrt(mean((v[m] .- o[m]) .^ 2)) : NaN
dstat(v, m) = any(m) ? (mean(v[m]), minimum(v[m]), maximum(v[m])) : (NaN, NaN, NaN)

# ── summary table ────────────────────────────────────────────────────────────
sets = [("TRAIN", TRAIN_DAYS), ("VAL", VAL_DAYS), ("TEST", TEST_DAYS), ("COMMON TEST", COMMON_TEST), ("ALL", collect(0:114))]
summ = Any[["config", "set", "T_hourly_RMSE", "T_5min_RMSE", "T_bias", "T_night_bias", "T_day_bias",
            "VP_day_bias_hPa", "VP_night_bias_hPa", "U8_bias"]]
println("\n", "="^118)
println("  FULL SEASON (days 0-114). T in K, VP in hPa. Measured air T is the reference.")
@printf("  %-5s %-12s | %7s %7s %7s %7s %7s | %7s %7s | %7s\n", "cfg", "set", "T hRMS", "T 5m", "bias", "night", "day",
        "VPday", "VPnite", "U8 bias")
for k in KEYS
    println("  ", "-"^116)
    T = R[k].T; VPm = VP(T, R[k].RH)
    for (sn, ds) in sets
        ins = in.(cd_, Ref(Set(ds))); m = okT .& ins; mv_ = okV .& ins; mu = okU .& ins
        row = [hourly_rmse(T, Tobs_g, m), rmse_(T, Tobs_g, m), mb(T, Tobs_g, m), mb(T, Tobs_g, m .& nite), mb(T, Tobs_g, m .& dayt),
               mb(VPm, VPo, mv_ .& dayt), mb(VPm, VPo, mv_ .& nite), mb(R[k].U, U8obs, mu)]
        @printf("  %-5s %-12s | %7.2f %7.2f %+7.2f %+7.2f %+7.2f | %+7.2f %+7.2f | %+7.3f\n", k, sn, row...)
        push!(summ, vcat(k, sn, round.(row, digits = 3)))
    end
end
writedlm(joinpath(FIGDIR, "AiCU_season_summary.csv"), permutedims(hcat(summ...)), ',')

# ── daily table ──────────────────────────────────────────────────────────────
hdr = ["day", "set", "meas_Tmean", "meas_Tmin", "meas_Tmax", "meas_VP", "meas_U8"]
for k in KEYS; append!(hdr, ["$(k)_Tmean", "$(k)_Tmin", "$(k)_Tmax", "$(k)_night_bias", "$(k)_day_bias", "$(k)_VP", "$(k)_U8"]); end
daily = Any[hdr]
for d in 0:114
    m = okT .& (cd_ .== d); mv_ = okV .& (cd_ .== d); mu = okU .& (cd_ .== d)
    mo = dstat(Tobs_g .- 273.15, m)
    row = Any[d, setname(d), mo..., any(mv_) ? mean(VPo[mv_]) : NaN, any(mu) ? mean(U8obs[mu]) : NaN]
    for k in KEYS
        T = R[k].T; append!(row, [dstat(T .- 273.15, m)..., mb(T, Tobs_g, m .& nite), mb(T, Tobs_g, m .& dayt),
                                  any(mv_) ? mean(VP(T, R[k].RH)[mv_]) : NaN, any(mu) ? mean(R[k].U[mu]) : NaN])
    end
    push!(daily, [x isa AbstractFloat ? round(x, digits = 3) : x for x in row])
end
writedlm(joinpath(FIGDIR, "AiCU_season_daily.csv"), permutedims(hcat(daily...)), ',')
D = permutedims(hcat(daily[2:end]...)); col(name) = Float64.(D[:, findfirst(==(name), hdr)]); dd = col("day")

# ── figures ──────────────────────────────────────────────────────────────────
CL = Dict("OL33" => :darkorange, "CL36" => :steelblue, "CL38" => :seagreen, "NV39" => :purple)
function shade!(pl)
    for (a, b) in TRAIN_BLOCKS; vspan!(pl, [a - 0.5, b + 0.5], color = :green, alpha = 0.08, label = ""); end
    for (a, b) in VAL_BLOCKS;   vspan!(pl, [a - 0.5, b + 0.5], color = :blue,  alpha = 0.08, label = ""); end
    pl
end
SHOWK = ["OL33", "CL38", "NV39"]
pls = Any[]
for (st, ttl) in (("Tmax", "daily maximum air T"), ("Tmean", "daily mean air T"), ("Tmin", "daily minimum air T"))
    pl = shade!(plot(ylabel = "°C", title = ttl * "  (green = TRAIN, blue = VAL, white = TEST)", titlefontsize = 9,
                     legend = st == "Tmax" ? :bottomleft : false, legendfontsize = 7))
    plot!(pl, dd, col("meas_" * st), lc = :black, lw = 2.5, label = "measured")
    for k in SHOWK; plot!(pl, dd, col("$(k)_" * st), lc = CL[k], lw = 1.4, label = LABS[k]); end
    push!(pls, pl)
end
savefig(plot(pls..., layout = (3, 1), size = (1200, 900), xlabel = "day of season"), joinpath(FIGDIR, "AiCU_season_T_daily.png"))

pb = Any[]
for (st, ttl) in (("night_bias", "night bias (model − measured), K"), ("day_bias", "daytime bias (model − measured), K"))
    pl = shade!(plot(ylabel = "K", title = ttl, titlefontsize = 9, legend = st == "night_bias" ? :topleft : false, legendfontsize = 7))
    hline!(pl, [0], lc = :gray, ls = :dash, label = "")
    for k in KEYS; plot!(pl, dd, col("$(k)_" * st), lc = CL[k], lw = 1.3, label = LABS[k]); end
    push!(pb, pl)
end
savefig(plot(pb..., layout = (2, 1), size = (1200, 650), xlabel = "day of season"), joinpath(FIGDIR, "AiCU_season_bias.png"))

pv = shade!(plot(ylabel = "hPa", title = "daily mean vapour pressure", titlefontsize = 9, legend = :topright, legendfontsize = 7))
plot!(pv, dd, col("meas_VP"), lc = :black, lw = 2.5, label = "measured")
for k in SHOWK; plot!(pv, dd, col("$(k)_VP"), lc = CL[k], lw = 1.3, label = LABS[k]); end
pu = shade!(plot(ylabel = "opening 0-1", title = "daily mean roof-vent opening (OL33 replays the log)", titlefontsize = 9, legend = false))
plot!(pu, dd, col("meas_U8"), lc = :black, lw = 2.5)
for k in ["CL38", "NV39"]; plot!(pu, dd, col("$(k)_U8"), lc = CL[k], lw = 1.3); end
savefig(plot(pv, pu, layout = (2, 1), size = (1200, 650), xlabel = "day of season"), joinpath(FIGDIR, "AiCU_season_vp_vent.png"))

SHOWD = parse.(Int, split(get(ENV, "SHOW_DAYS", "10,24,45,56,73,93,105,112"), ","))
pd = Any[]
for d in SHOWD
    c = [col_of(d, j) for j in 0:NPD-1]; h = (0:NPD-1) .* (300 / 3600)
    pl = plot(h, Tobs_g[c] .- 273.15, lc = :black, lw = 2.5, label = "measured", title = "day $d ($(setname(d)))",
              titlefontsize = 9, legend = d == SHOWD[1] ? :topleft : false, legendfontsize = 6, xlabel = "hour", ylabel = "°C")
    for k in SHOWK; plot!(pl, h, R[k].T[c] .- 273.15, lc = CL[k], lw = 1.2, label = LABS[k]); end
    push!(pd, pl)
end
savefig(plot(pd..., layout = (2, 4), size = (1600, 700)), joinpath(FIGDIR, "AiCU_season_days.png"))
println("\nFigures and CSVs written to ", FIGDIR, ":")
for f in ("AiCU_season_T_daily.png", "AiCU_season_bias.png", "AiCU_season_vp_vent.png", "AiCU_season_days.png",
          "AiCU_season_daily.csv", "AiCU_season_summary.csv"); println("  ", f); end
