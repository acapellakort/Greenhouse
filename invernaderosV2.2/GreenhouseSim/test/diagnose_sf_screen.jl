# =============================================================================
# diagnose_sf_screen.jl -- can physical parameter values close the biases of the SF structure? (Note 46 -> 47)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1-2 min, no MCMC):
#     julia --project=. test/diagnose_sf_screen.jl
#
# Structure SF of Note 46: Note 44 D3 (C3 flags, gamma1 x 2, soil_dyn 14 d) + floor_full + solar_split + deep_k.
# Starting point = Note 38 posterior. One parameter at a time:
#   gamma4   stomatal resistance multiplier x1, 1.5, 2, 3, 5   (x1 = Note 38, at its lower bound = max transpiration)
#   cov_abs  0.05, 0.10, 0.20                                   (sun absorbed by the cover)
#   k_deep   1, 2, 5, 10 W m-2 K-1                               (floor -> deep soil)
#   C_soil   1.5e5, 3e5, 6e5, 1.2e6 J m-2 K-1                    (active floor layer ~0.1 / 0.2 / 0.4 / 0.8 m)
#   gamma1   pipe length x2 (D3), x3, x4 of the JSON value
#   epsil6   0.44 (V2.2), 0.22, 0                               (air -> cover radiation lump, Note 32)
#   sd_tau   7, 14, 30 d
# Then two combined runs built from the sweeps (printed at the end). k_ground stays at the Note 38 value.
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
function row(lab, p)
    s = scores(p)
    @printf("  %-28s | %6.3f | %+6.2f %+6.2f %+7.2f %+7.2f %+7.2f | %+6.2f %+6.2f %+6.2f %+7.3f | %5.2f\n",
            lab, s.h, s.nb..., s.d, s.vpd, s.vpn, s.ud, s.amp)
    s
end
println("\n", "="^124)
@printf("  %-28s | %6s | %6s %6s %7s %7s %7s | %6s %6s %6s %7s | %5s\n", "run", "hRMS", "n0-8", "n9-48", "n49-100",
        "n101-9", "n110-4", "T day", "VPday", "VPnit", "U8 day", "T5amp")
println("  ", "-"^122)
ref = row("R0 D3 reference (old)", pD3)
row("SF0 Note 38 params", SF0)
up1(p, k, v) = update_params(p, [k], [v])
println("  -- gamma4 (x Note 38) --")
for m in (1.5, 2.0, 3.0, 5.0); row("gamma4 x$(m)", up1(SF0, "gamma4", m * SF0.gamma4)); end
println("  -- cov_abs --")
for v in (0.05, 0.20); row("cov_abs $(v)", up1(SF0, "cov_abs", v)); end
println("  -- k_deep --")
for v in (1.0, 5.0, 10.0); row("k_deep $(v)", up1(SF0, "k_deep", v)); end
println("  -- C_soil --")
for v in (3e5, 6e5, 1.2e6); row("C_soil $(v)", up1(SF0, "C_soil", v)); end
println("  -- gamma1 (x JSON) --")
for m in (3.0, 4.0); row("gamma1 x$(m)", up1(SF0, "gamma1", m * params.gamma1)); end
println("  -- epsil6 --")
for v in (0.22, 0.0); row("epsil6 $(v)", up1(SF0, "epsil6", v)); end
println("  -- sd_tau --")
for v in (7.0, 30.0); row("sd_tau $(v) d", up1(SF0, "sd_tau", v * 86400.0)); end
println("  -- combinations --")
row("C1 gamma4x3 C_soil6e5",            update_params(SF0, ["gamma4", "C_soil"], [3 * SF0.gamma4, 6e5]))
row("C2 C1 + gamma1x3",                 update_params(SF0, ["gamma4", "C_soil", "gamma1"], [3 * SF0.gamma4, 6e5, 3 * params.gamma1]))
row("C3 C2 + k_deep 5",                 update_params(SF0, ["gamma4", "C_soil", "gamma1", "k_deep"], [3 * SF0.gamma4, 6e5, 3 * params.gamma1, 5.0]))
println("""
  How to read it:
   - Look for parameters that move the DAY bias up (SF0: about -1.4 K) without warming VP day further, and
     parameters that lift the NIGHTS (SF0: -0.5 to -2.7 K).
   - A sweep that has no effect means that parameter cannot be the missing piece.
   - If physical values (C_soil up to ~1e6, gamma1 up to the real pipe length, k_deep 1-10) cannot bring the
     combination near 1.3-1.5 K, a heat input is still missing and the calibration would just find new
     compensations; then we look for that input before the long run.""")
