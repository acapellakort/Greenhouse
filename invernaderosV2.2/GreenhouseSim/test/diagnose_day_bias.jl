# =============================================================================
# diagnose_day_bias.jl -- where does the remaining ~1 K daytime cold bias of the SF structure live? (Note 47)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, < 1 min, no MCMC):
#     julia --project=. test/diagnose_day_bias.jl
#
# Two runs: R0 (Note 44 D3, old structure) and Q (best physical candidate of screen 2:
# SF + gamma4 = 82 + nu1 = 0.65 + C_soil = 6e5 + k_ground = 5). Daytime 5-min samples (Iglob > 50 W m-2) are
# binned by outside radiation, logged vent opening, energy screen, hour of day and season period.
# Reading:
#   bias grows with radiation, vents mostly open      -> solar input or a sun-heated sensor
#   bias concentrated where logged vents are open     -> vent controller offset/gain (calibration can fix it)
#   bias concentrated with vents closed, low sun      -> a missing heat input (pipes, lamps) or too-strong losses
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
TR = Dict{String,Any}()
for (nm_, p) in (("R0", pD3), ("Q", Q))
    us_ = simrun(p, weather_c, ctrl_v, u0_9)
    TR[nm_] = (T = us_[2, :], U = clamp.(us_[7, :], 0.0, 1.0))
end
Igc  = A(Ig_v).(predtimes)
U1c  = A(U1v).(predtimes)
hrc  = [floor(Int, ((c - 2) % NPD) * 300 / 3600) for c in 1:NCOL]; hrc[1] = 0
D = okT .& dayt
function tab(title, labels, masks)
    println("\n  ", title)
    @printf("  %-22s | %6s | %8s %8s | %8s %8s | %8s\n", "bin", "n", "R0 bias", "Q bias", "R0 dU8", "Q dU8", "logged U8")
    for (lab, m) in zip(labels, masks)
        mm = D .& m; mu = mm .& okU
        n = count(mm); n < 50 && continue
        @printf("  %-22s | %6d | %+8.2f %+8.2f | %+8.3f %+8.3f | %8.3f\n", lab, n,
                mb(TR["R0"].T, Tobs_g, mm), mb(TR["Q"].T, Tobs_g, mm),
                mb(TR["R0"].U, U8obs, mu), mb(TR["Q"].U, U8obs, mu), any(mu) ? mean(U8obs[mu]) : NaN)
    end
end
tab("by outside global radiation (W m-2)", ["50-150", "150-300", "300-500", "500-700", ">700"],
    [(Igc .>= a) .& (Igc .< b) for (a, b) in ((50, 150), (150, 300), (300, 500), (500, 700), (700, 1e9))])
U8f = replace(U8obs, NaN => -1.0)
tab("by LOGGED vent opening", ["closed (<0.02)", "0.02-0.10", "0.10-0.30", ">0.30"],
    [(U8f .>= 0) .& (U8f .< 0.02), (U8f .>= 0.02) .& (U8f .< 0.10), (U8f .>= 0.10) .& (U8f .< 0.30), U8f .>= 0.30])
tab("by energy screen", ["open (U1<0.1)", "partly", "closed (U1>0.9)"],
    [U1c .< 0.1, (U1c .>= 0.1) .& (U1c .<= 0.9), U1c .> 0.9])
tab("by hour of day", ["06-09", "09-12", "12-15", "15-18", "18-21"],
    [(hrc .>= a) .& (hrc .< b) for (a, b) in ((6, 9), (9, 12), (12, 15), (15, 18), (18, 21))])
tab("by season period (days)", ["0-8", "9-48", "49-100", "101-109", "110-114"],
    [inset(collect(a:b)) for (a, b) in PER])
@printf("\n  all daytime: R0 %+.2f K, Q %+.2f K;  Q hourly RMSE %.3f K\n", mb(TR["R0"].T, Tobs_g, D), mb(TR["Q"].T, Tobs_g, D),
        hourly_rmse(TR["Q"].T, Tobs_g, okT .& ALLD))
println("""
  dU8 = model minus logged vent opening (hourly means not used here: 5-min samples).
  Note: hour of day assumes the season starts at midnight (as the hourly bins elsewhere).""")
