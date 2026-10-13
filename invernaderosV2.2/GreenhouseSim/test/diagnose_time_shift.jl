# =============================================================================
# diagnose_time_shift.jl -- is the morning-warm / afternoon-cold daytime pattern a clock offset? (Note 48)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     julia --project=. test/diagnose_time_shift.jl
#
# Note 47 day-bias table: in BOTH structures the model is too warm in the morning and too cold in the afternoon
# (R0 +1.47 K at 06-09 -> -0.46 K at 15-18; Q +0.20 -> -1.77 K). Two explanations:
#   (a) a clock offset between the weather inputs and the greenhouse logs (UTC vs local time, summer time),
#   (b) missing thermal inertia (the real greenhouse stores heat in the morning and releases it later).
# Part A (data only): timing of the outside-radiation centroid in the log clock, against solar noon at
#   Bleiswijk (4.5 E): ~11.7 h in UTC, ~12.7 h in CET (UTC+1), ~13.7 h in CEST (UTC+2). Also the timing of the
#   measured inside-air maximum, and whether the season starts at midnight.
# Part B (forward): structure Q of Note 47 with all weather inputs delayed by -2, -1, 0, +1, +2 h.
#   If one shift flattens the hour-of-day bias, (a); if none does, (b).
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n37"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
for k in (:hc_log_pipe, :vc_min_log, :soil_dyn, :sd_tau, :floor_full, :solar_split, :cov_abs, :deep_k, :k_deep)
    hasproperty(params, k) || error("constants_climate.json has no $k entry (Note 42-46 flags missing)")
end
isdefined(GreenhouseSim, :floor_full_fluxes) || error("src/climate.jl has no floor_full_fluxes (Note 46 edit missing)")

parts = Vector{Vector{Float64}}()
for s in (1, 2)
    d = JLD2.load(joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n37_seed$(s).jld2"))
    push!(parts, vec(mean(d["chain"][d["burn_in"]+1:end, 1:NM], dims = 1)))
end
pm = (parts[1] .+ parts[2]) ./ 2
p38 = update_params(params_cl, CQ_MODEL, pm)

VentMinv = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    vm = filter(!isnan, [_num(vip[i, Symbol("VentLeewMin_Vip_$z")]) for z in 1:6])
    push!(VentMinv, isempty(vm) ? 0.0 : clamp(mean(vm) / 200, 0.0, 1.0))
end
length(VentMinv) == length(ts) || error("row filter mismatch (VentMin)")
ctrl_v = ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                        Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv))
day0 = [i for i in eachindex(ts) if ts[i] < t0 + 86400 && !isnan(Tair_obs[i])]
u0_9 = vcat(u0_8, mean(Tair_obs[day0]))

function simrun(p, w, ctrl, u0x)
    us = Matrix{Float64}(undef, length(u0x), 1 + length(days) * NPD); us[:, 1] = u0x
    u = copy(u0x); colx = 2
    for d in days
        ctx = GreenhouseSim.SimContext(p, w, ctrl, t0 + 86400.0 * d, 86400.0 * d)
        sol = GreenhouseSim.solve(GreenhouseSim.ODEProblem(GreenhouseSim.rhs!, u, tspan, ctx), GreenhouseSim.Tsit5();
                                  saveat = teval, reltol = RELTOL)
        length(sol.t) == NPD || error("solver stopped on day $d")
        S = Array(sol); us[:, colx:colx+NPD-1] = S; colx += NPD; u = S[:, end]
    end
    us
end
pD3 = update_params(p38, ["hc_log_pipe", "vc_min_log", "gamma1", "soil_dyn", "sd_tau"],
                    [1.0, 1.0, 2.0 * p38.gamma1, 1.0, 14 * 86400.0])
fl(p; ff = 0.0, ss = 0.0, dk = 0.0, kg = p.k_ground) =
    update_params(p, ["floor_full", "solar_split", "deep_k", "k_ground"], [ff, ss, dk, kg])
@printf("\nNote 38 posterior: k_ground = %.1f, vc_gain = %.2f; fixed: k_deep = %.1f, cov_abs = %.2f, C_soil = %.2g\n",
        p38.k_ground, p38.vc_gain, p38.k_deep, p38.cov_abs, p38.C_soil)

okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g); okU = .!isnan.(U8obs); okU[1] = false
inset(ds) = in.(cd_, Ref(Set(ds)))
ALLD = inset(collect(0:114)); COM = inset(COMMON_TEST)
VP(T, RH) = RH ./ 100 .* Pws.(T) ./ 100; VPo = VP(Tobs_g, RHobs_g)
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
function hourly_rmse(v, o, mask)
    e = Float64[]
    for cols in values(HB); k = [c for c in cols if mask[c]]; length(k) >= 6 && push!(e, mean(v[k]) - mean(o[k])); end
    isempty(e) ? NaN : sqrt(mean(abs2, e))
end
PER = [(0, 8), (9, 48), (49, 100), (101, 109), (110, 114)]
famp(T5) = mean([maximum(T5[cd_ .== d]) - minimum(T5[cd_ .== d]) for d in 1:114])


SF0 = update_params(pD3, ["floor_full", "solar_split", "deep_k"], [1.0, 1.0, 1.0])
@printf("\nStarting values: gamma4 = %.3g, cov_abs = %.2f, k_deep = %.1f, C_soil = %.2g, gamma1 = %.2f (JSON %.2f), epsil6 = %.2f, k_ground = %.1f\n",
        SF0.gamma4, SF0.cov_abs, SF0.k_deep, SF0.C_soil, SF0.gamma1, params.gamma1, SF0.epsil6, SF0.k_ground)

function scores(p)
    us_ = simrun(p, weather_c, ctrl_v, u0_9); T = us_[2, :]; RH = rhf.(T, us_[4, :]); U = clamp.(us_[7, :], 0.0, 1.0)
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    (h = hourly_rmse(T, Tobs_g, okT .& ALLD), nb = nb, d = mb(T, Tobs_g, okT .& ALLD .& dayt),
     vpd = mb(VP(T, RH), VPo, okV .& ALLD .& dayt), vpn = mb(VP(T, RH), VPo, okV .& ALLD .& nite),
     ud = mb(U, U8obs, okU .& ALLD .& dayt), amp = famp(us_[6, :]))
end

Q = update_params(SF0, ["gamma4", "nu1", "C_soil", "k_ground"], [82.0, 0.65, 6e5, 5.0])

# ── Part A: clocks ───────────────────────────────────────────────────────────
@printf("\nt0 = %.1f s;  t0/86400 = %.4f  (fraction of day at the season start: %.3f -> %s)\n",
        t0, t0 / 86400, mod(t0 / 86400, 1.0), mod(t0 / 86400, 1.0) < 1e-3 ? "starts at midnight" : "NOT midnight: hour bins are shifted")
dayof(t) = floor(Int, (t - t0) / 86400)
hod(t)   = mod(t - t0, 86400.0) / 3600 + 24 * mod(t0 / 86400, 1.0)
cent = Float64[]; tmax = Float64[]; dlist = Int[]
for d in 0:114
    k = [i for i in eachindex(ts) if dayof(ts[i]) == d]
    isempty(k) && continue
    I = Ig_v[k]; h = mod.(hod.(ts[k]), 24.0)
    sum(I) > 0 || continue
    push!(cent, sum(I .* h) / sum(I)); push!(dlist, d)
    Tk = Tair_obs[k]; ok = .!isnan.(Tk)
    push!(tmax, any(ok) ? h[ok][argmax(Tk[ok])] : NaN)
end
println("\n  PART A  radiation-weighted hour (centroid) and hour of the measured inside-T maximum, in the log clock")
@printf("  %-10s | %8s | %10s | %s\n", "days", "I centr.", "Tin max", "expected I centroid: UTC ~11.7, CET ~12.7, CEST ~13.7")
for (a, b) in ((0, 20), (21, 40), (41, 60), (61, 80), (81, 100), (101, 114))
    m = (dlist .>= a) .& (dlist .<= b)
    any(m) || continue
    @printf("  %3d-%-6d | %8.2f | %10.2f |\n", a, b, mean(cent[m]), mean(filter(!isnan, tmax[m])))
end
println("  (a jump of ~1 h in the centroid between two periods = a summer/winter time change inside the log)")

# ── Part B: delayed weather ──────────────────────────────────────────────────
As(v, s) = stepfun(ts .+ s, v)
wshift(s) = WeatherInputs(As(Ig_v, s), As(Tsky_v, s), As(Tout_v, s), weather_c.TmechCool, weather_c.Tsoil,
                          As(Wind_v, s), weather_c.CO2out, As(VPout_v, s), As(Ig_v, s), weather_c.LAI)
Q = update_params(SF0, ["gamma4", "nu1", "C_soil", "k_ground"], [82.0, 0.65, 6e5, 5.0])
hrc = [floor(Int, mod(((c - 2) % NPD) * 300 / 3600 + 24 * mod(t0 / 86400, 1.0), 24.0)) for c in 1:NCOL]; hrc[1] = 0
HB_ = ((6, 9), (9, 12), (12, 15), (15, 18), (18, 21))
println("\n  PART B  structure Q, weather delayed by s hours. Day bias (K) by hour of day; hRMS all days; night bias 49-100.")
@printf("  %-8s | %6s | %s | %7s | %7s %7s\n", "shift", "hRMS", join([@sprintf("%7s", "$(a)-$(b)") for (a, b) in HB_], " "),
        "n49-100", "T day", "U8 day")
for s in (-2.0, -1.0, 0.0, 1.0, 2.0)
    us_ = simrun(Q, s == 0.0 ? weather_c : wshift(s * 3600.0), ctrl_v, u0_9)
    T = us_[2, :]; U = clamp.(us_[7, :], 0.0, 1.0)
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    @printf("  %+5.0f h  | %6.3f | %s | %+7.2f | %+7.2f %+7.3f\n", s, hourly_rmse(T, Tobs_g, okT .& ALLD),
            join([@sprintf("%+7.2f", x) for x in hb], " "), mb(T, Tobs_g, okT .& nite .& inset(collect(49:100))),
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(U, U8obs, okU .& ALLD .& dayt))
end
println("""
  How to read it:
   - Part A: the centroid should sit at ONE of 11.7 / 12.7 / 13.7 h. A jump of 1 h within the season means the log
     switches between summer and winter time.
   - Part B: the shift whose hour-of-day row is flattest is the clock offset between weather and greenhouse logs.
     The 0 h row must reproduce Note 47 Q (hRMS 1.698). If no shift flattens the row, the pattern is missing
     thermal inertia, not timing.
   - Note: the measured inside T and the controls (vents, pipe, screen) are in the greenhouse-log clock and are
     not shifted; only the outside weather is.""")
