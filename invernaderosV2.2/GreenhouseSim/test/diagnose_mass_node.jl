# =============================================================================
# diagnose_mass_node.jl -- does a lumped internal thermal mass remove the morning-warm / afternoon-cold drift? (Note 48)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     julia --project=. test/diagnose_mass_node.jl
#
# The time-shift test ruled out a clock offset (the log is in local time, CEST -> CET at day 75 = 28 Oct 2018,
# weather and greenhouse in the same clock). Remaining explanation for the ~2 K morning->afternoon drift: missing
# thermal inertia. New default-off flag mass_node (src/forward_map.jl, backup *.bak_premass):
#   C_mass dTm/dt = h_mass (T2 - Tm) + f_mass * Iglob,    air loses h_mass (T2 - Tm);  state u[10]
# (structure, substrate slabs, crop water, pipe water; ~5e4 J m-2 K-1 = 12 L water per m2).
# Runs on structure Q of Note 47:
#   Q            no mass (9 states);  check: 10 states, flag off (must equal Q)
#   M C/h        C_mass 2e4 / 5e4 / 1e5 / 2e5  x  h_mass 3 / 10, no solar share
#   MS           C_mass 5e4 or 1e5, h_mass 5, f_mass 0.1 / 0.2 (the mass also absorbs sun: lagged solar heat)
#   instant sun  tau1 0.85 (more transmitted sun), eta1 0.2 (more sun absorbed straight into the air) -- no lag
# A good candidate flattens the hour-of-day row AND lifts the day without warming the nights.
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


Q   = update_params(SF0, ["gamma4", "nu1", "C_soil", "k_ground"], [82.0, 0.65, 6e5, 5.0])
u0_10 = vcat(u0_9, u0_9[2])
hrc = [floor(Int, mod(((c - 2) % NPD) * 300 / 3600, 24.0)) for c in 1:NCOL]; hrc[1] = 0
HB_ = ((6, 9), (9, 12), (12, 15), (15, 18), (18, 21))
function mrow(lab, p, u0x)
    us_ = simrun(p, weather_c, ctrl_v, u0x); T = us_[2, :]; RH = rhf.(T, us_[4, :]); U = clamp.(us_[7, :], 0.0, 1.0)
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    spread = maximum(hb) - minimum(hb)
    @printf("  %-26s | %6.3f | %s | %4.2f | %s | %+6.2f %+6.2f %+7.3f\n", lab, hourly_rmse(T, Tobs_g, okT .& ALLD),
            join([@sprintf("%+6.2f", x) for x in hb], " "), spread, join([@sprintf("%+6.2f", x) for x in nb], " "),
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& dayt), mb(U, U8obs, okU .& ALLD .& dayt))
    T
end
println("\n", "="^150)
@printf("  %-26s | %6s | %s | %4s | %s | %6s %6s %7s\n", "run", "hRMS", join([@sprintf("%6s", "$(a)-$(b)") for (a, b) in HB_], " "),
        "sprd", join([@sprintf("%6s", "n$(a)") for (a, b) in PER], " "), "T day", "VPday", "U8 day")
println("  ", "-"^148)
Tq  = mrow("Q (9 states)", Q, u0_9)
Tq0 = mrow("check: 10 states, mass off", Q, u0_10)
@printf("  [check 0] max |dT| between the two rows: %.2e K\n", maximum(abs.(Tq .- Tq0)))
M(p; C = 5e4, h = 5.0, f = 0.0) = update_params(p, ["mass_node", "C_mass", "h_mass", "f_mass"], [1.0, C, h, f])
println("  -- mass, no solar share --")
for C in (2e4, 5e4, 1e5, 2e5), h in (3.0, 10.0)
    mrow("M C=$(round(Int, C/1000))k h=$(round(Int, h))", M(Q; C = C, h = h), u0_10)
end
println("  -- mass that also absorbs sun --")
for C in (5e4, 1e5), f in (0.1, 0.2)
    mrow("MS C=$(round(Int, C/1000))k h=5 f=$(f)", M(Q; C = C, h = 5.0, f = f), u0_10)
end
println("  -- more sun, no lag --")
mrow("tau1 0.85", update_params(Q, ["tau1"], [0.85]), u0_9)
mrow("eta1 0.2", update_params(Q, ["eta1"], [0.2]), u0_9)
println("""
  Columns: day bias (K) by hour of day, sprd = max - min of those five (Q: ~2.0 K), night bias by period,
  day T bias, VP day (hPa), U8 day (model - logged).
  How to read it:
   - A mass that lowers sprd clearly (towards < 1 K) is the missing inertia; f_mass > 0 then also lifts the day.
   - If only the 'more sun' rows lift the day and sprd stays ~2 K, the drift is something else (e.g. afternoon
     transpiration or vent dynamics), not inertia.""")
