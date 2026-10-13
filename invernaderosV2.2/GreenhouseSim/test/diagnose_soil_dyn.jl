# =============================================================================
# diagnose_soil_dyn.jl -- replace the measured 7-day soil anchor by a model deep-soil state (Note 43)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1-2 min, no MCMC):
#     julia --project=. test/diagnose_soil_dyn.jl
#
# All runs: Note 38 posterior mean, closed-loop vents, C3 structure of Note 43
#   (hc_log_pipe = 1, vc_min_log = 1 with the logged VentLeewMin input, pipe length gamma1 x 2).
#   R0   current model (Note 38), for reference
#   D0   C3 with the MEASURED 7-day soil anchor (Note 43 reference; uses observations as input)
#   D1   C3 with soil_dyn, sd_tau = 3 days
#   D2   C3 with soil_dyn, sd_tau = 7 days
#   D3   C3 with soil_dyn, sd_tau = 14 days
#   D4   D2 with k_ground = 10 (is the bound-pinned coupling still needed?)
#   D5   D2 with the deep soil starting at 18 °C instead of the day-0 mean air T
# soil_dyn: dTdeep/dt = (T_floor - Tdeep)/sd_tau, state u[9]; the floor couples to Tdeep through k_ground.
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n37"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
for k in (:hc_log_pipe, :vc_min_log, :soil_dyn, :sd_tau)
    hasproperty(params, k) || error("constants_climate.json has no $k entry (Note 42/43 flags missing)")
end

parts = Vector{Vector{Float64}}()
for s in (1, 2)
    d = JLD2.load(joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n37_seed$(s).jld2"))
    push!(parts, vec(mean(d["chain"][d["burn_in"]+1:end, 1:NM], dims = 1)))
end
pm = (parts[1] .+ parts[2]) ./ 2
p38 = update_params(params_cl, CQ_MODEL, pm)

# inputs (as diagnose_night_heating.jl)
VentMinv = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    vm = filter(!isnan, [_num(vip[i, Symbol("VentLeewMin_Vip_$z")]) for z in 1:6])
    push!(VentMinv, isempty(vm) ? 0.0 : clamp(mean(vm) / 200, 0.0, 1.0))
end
length(VentMinv) == length(ts) || error("row filter mismatch (VentMin)")
Tfill = copy(Tair_obs); for i in eachindex(Tfill); isnan(Tfill[i]) && (Tfill[i] = i > 1 ? Tfill[i-1] : 291.15); end
soil7 = similar(Tfill); j0 = 1; acc = 0.0
for i in eachindex(ts)
    global acc += Tfill[i]
    while ts[j0] < ts[i] - 7 * 86400; global acc -= Tfill[j0]; global j0 += 1; end
    soil7[i] = acc / (i - j0 + 1)
end
weather_s = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, A(soil7),
                          A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), A(LAI_row))
ctrl_v = ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8z), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                        Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv))

day0 = [i for i in eachindex(ts) if ts[i] < t0 + 86400 && !isnan(Tair_obs[i])]
Tdeep0 = mean(Tair_obs[day0])
@printf("\nDeep-soil start: mean measured air T on day 0 = %.2f °C (D5 uses 18.0 °C)\n", Tdeep0 - 273.15)
u0_9  = vcat(u0_8, Tdeep0)
u0_9c = vcat(u0_8, 18.0 + 273.15)

function simrun(p, w, ctrl, u0x)
    us = Matrix{Float64}(undef, length(u0x), 1 + length(days) * NPD); us[:, 1] = u0x
    u = copy(u0x); col = 2
    for d in days
        ctx = GreenhouseSim.SimContext(p, w, ctrl, t0 + 86400.0 * d, 86400.0 * d)
        sol = GreenhouseSim.solve(GreenhouseSim.ODEProblem(GreenhouseSim.rhs!, u, tspan, ctx), GreenhouseSim.Tsit5();
                                  saveat = teval, reltol = RELTOL)
        length(sol.t) == NPD || error("solver stopped on day $d")
        S = Array(sol); us[:, col:col+NPD-1] = S; col += NPD; u = S[:, end]
    end
    us
end
pC3 = update_params(p38, ["hc_log_pipe", "vc_min_log", "gamma1"], [1.0, 1.0, 2.0 * p38.gamma1])
pSD(tau; kg = p38.k_ground) = update_params(pC3, ["soil_dyn", "sd_tau", "k_ground"], [1.0, tau * 86400.0, kg])

# check 0: 9 states with soil_dyn off vs 8 states (only the adaptive step norm differs)
ua = simrun(p38, weather_c, ctrl_x, u0_8); ub = simrun(p38, weather_c, ctrl_x, u0_9)
@printf("[check 0] 9 states, soil_dyn off vs 8 states: max |ΔT| = %.3e K, mean |ΔT| = %.3e K (step-size effect only)\n",
        maximum(abs.(ua[2, :] .- ub[2, :])), mean(abs.(ua[2, :] .- ub[2, :])))

runs = [("R0 current (N38)",       p38,         weather_c, ctrl_x, u0_8),
        ("D0 C3, measured anchor", pC3,         weather_s, ctrl_v, u0_8),
        ("D1 soil_dyn 3 d",        pSD(3),      weather_c, ctrl_v, u0_9),
        ("D2 soil_dyn 7 d",        pSD(7),      weather_c, ctrl_v, u0_9),
        ("D3 soil_dyn 14 d",       pSD(14),     weather_c, ctrl_v, u0_9),
        ("D4 7 d, k_ground 10",    pSD(7; kg = 10.0), weather_c, ctrl_v, u0_9),
        ("D5 7 d, start 18 °C",    pSD(7),      weather_c, ctrl_v, u0_9c)]
res = Any[]
println("Running $(length(runs)) full-season simulations ...")
for (nm_, p, w, ct, u0x) in runs
    tt = @elapsed (us_ = simrun(p, w, ct, u0x))
    push!(res, (nm_, us_[2, :], rhf.(us_[2, :], us_[4, :]), clamp.(us_[7, :], 0.0, 1.0),
                length(u0x) >= 9 ? us_[9, :] : fill(NaN, size(us_, 2))))
    @printf("  %-24s %.1f s\n", nm_, tt)
end

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

println("\n", "="^150)
println("  FULL SEASON. hRMS = hourly T RMSE (K). Night/day bias (K) by period; VP night (hPa) all days. Tdeep = deep-soil state (°C) on days 0/30/60/90/114.")
@printf("  %-24s | %6s %6s | %7s %7s %8s %7s %7s | %6s %7s | %s\n", "run", "hRMS", "hRMS", "n 0-8", "n 9-48", "n 49-100",
        "n101-9", "n110-4", "day", "VPnite", "Tdeep d0/30/60/90/114")
@printf("  %-24s | %6s %6s | %43s | %6s %7s |\n", "", "all", "common", "", "all", "all")
println("  ", "-"^148)
for (nm_, T, RH, U, Td) in res
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    tds = all(isnan, Td) ? "   (fixed anchor)" : join([@sprintf("%5.1f", Td[col_of(d, 144)] - 273.15) for d in (0, 30, 60, 90, 114)], " ")
    @printf("  %-24s | %6.3f %6.3f | %+7.2f %+7.2f %+8.2f %+7.2f %+7.2f | %+6.2f %+7.2f | %s\n", nm_,
            hourly_rmse(T, Tobs_g, okT .& ALLD), hourly_rmse(T, Tobs_g, okT .& COM), nb...,
            mb(T, Tobs_g, okT .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& nite), tds)
end
@printf("  measured 7-day anchor (D0 input), days 0/30/60/90/114: %s\n",
        join([@sprintf("%5.1f", A(soil7)(t0 + 86400.0 * d + 43200) - 273.15) for d in (0, 30, 60, 90, 114)], " "))
println("""
  How to read it:
   - D1-D3 close to D0           -> a model deep-soil state reproduces the measured anchor; no observation leak needed.
   - D4 close to D2              -> k_ground no longer needs to sit at its bound.
   - D5 vs D2                    -> sensitivity to the unknown initial soil temperature (spin-up).
   Parameters are the Note 38 posterior with gamma1 x 2; this is a sensitivity test, not a calibration.""")

pl = plot(ylabel = "night bias (K)", xlabel = "day of season", title = "Daily night bias, model − measured", legend = :bottomleft,
          legendfontsize = 7, size = (1200, 450))
for (a, b) in TRAIN_BLOCKS; vspan!(pl, [a - 0.5, b + 0.5], color = :green, alpha = 0.08, label = ""); end
hline!(pl, [0], lc = :gray, ls = :dash, label = "")
for (k, (nm_, T, _, _, _)) in enumerate(res)
    k in (1, 2, 4, 6, 7) || continue
    plot!(pl, 0:114, [mb(T, Tobs_g, okT .& nite .& (cd_ .== d)) for d in 0:114], lw = k == 1 ? 2.2 : 1.4, label = nm_)
end
fig = joinpath(FIGDIR, "AiCU_soil_dyn.png"); savefig(pl, fig); println("Figure: ", fig)
