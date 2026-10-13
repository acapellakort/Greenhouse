# =============================================================================
# diagnose_heating_bias.jl -- why does the Vanthoor cold bias grow over the season?
# =============================================================================
# Run from the GreenhouseSim folder (~1 min, no MCMC):
#     julia --project=. test/diagnose_heating_bias.jl
#
# Full-season (days 0-114) forward runs, no inference:
#   R1  B         joint posterior mean theta_v (Note 31), LAI = 2.0 (as calibrated)
#   R2  N29       nu4=4.46e-5, gamma4=44, k_ground=24.5, U8_min_leak=0.05, LAI = 2.0
#   R3  N29+LAI   R2 with Vanthoor LAI = derived observed LAI(t) instead of 2.0
#   R4  B+LAI     R1 with derived observed LAI(t)
#   R5  N29+pipe  R2 with +PIPE_TEST_K on rows where Tpipe is the heating ESTIMATE
# Splits the bias into night (Iglob < 1 W/m2) / day, by week, and by Tpipe source
# (measured pipe / heating estimate / cold fallback), and correlates the daily
# night bias with outdoor T and with how often the heating estimate was used.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = get(ENV, "N_FULL", "115")
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const PIPE_TEST_K = parse(Float64, get(ENV, "PIPE_TEST_K", "10"))
const DG_FIG = joinpath(FIGDIR, "AiCU_heating_bias")

# ── Tpipe source per data row (same rules and row filter as the joint script) ─
src = Int[]            # 1 = measured pipe, 2 = heating estimate, 3 = cold fallback
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    To = _num(mt.Tout[i]); To_val = isnan(To) ? 15.0 : To
    pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    hv = filter(!isnan, [_num(vip.HeatTemp_Vip_1[i]), _num(vip.HeatTemp_Vip_2[i]), _num(vip.HeatTemp_Vip_3[i]),
                         _num(vip.HeatTemp_Vip_5[i]), _num(vip.HeatTemp_Vip_6[i])])
    Th = isempty(hv) ? 20.0 : mean(hv)
    valid = (!isnan(pg) && pg > PIPE_MIN_VALID) || (!isnan(pl) && pl > PIPE_MIN_VALID)
    push!(src, valid ? 1 : (To_val < Th - 2.0 && To_val < PIPE_WARM_TOUT) ? 2 : 3)
end
length(src) == length(ts) || error("row filter mismatch: $(length(src)) vs $(length(ts))")

# ── alternative inputs ───────────────────────────────────────────────────────
# derived observed LAI(t): LAI0 = 0.5 at day 0, then CropManagement-derived points, flat after
lai_d = vcat(0.0, Float64.(cm_days_all)); lai_v = vcat(0.5, cm_L_all)
function lai_of_day(x)
    x <= lai_d[1] && return lai_v[1]; x >= lai_d[end] && return lai_v[end]
    i = searchsortedlast(lai_d, x)
    return lai_v[i] + (x - lai_d[i]) / (lai_d[i+1] - lai_d[i]) * (lai_v[i+1] - lai_v[i])
end
laifun = tcal -> lai_of_day((tcal - t0) / 86400.0)
weather_lai = WeatherInputs(A(Ig_v), A(Tsky_v), A(Tout_v), _t -> 273.15, _t -> 18.0 + 273.15,
                            A(Wind_v), _t -> 834.7, A(VPout_v), A(Ig_v), laifun)
base_lai  = SimBase(params, weather_lai, controls, t0)
Tpipe_alt = [src[i] == 2 ? Tpipev[i] + PIPE_TEST_K : Tpipev[i] for i in eachindex(Tpipev)]
controls_pipe = ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8v), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_alt))
base_pipe = SimBase(params, weather, controls_pipe, t0)

θN29 = [4.46e-5, 44.0, 24.5, 0.05]
jf = sort(filter(f -> startswith(f, "AiCU_joint_optionB_d42_seed") && endswith(f, ".jld2") && !occursin("ckpt", f), readdir(POSTDIR)))
θB = if isempty(jf)
    @warn "no joint chains found; R1/R4 use N29 values"; θN29
else
    P = vcat([let c = JLD2.load(joinpath(POSTDIR, f)); c["chain"][c["burn_in"]+1:end, IV] end for f in jf]...)
    vec(mean(P, dims = 1))
end
println("\ntheta_v (nu4, gamma4, k_ground, U8_min_leak):  B = ", round.(θB, sigdigits = 3), "   N29 = ", θN29)

runs = [("R1 B",        θB,   base),
        ("R2 N29",      θN29, base),
        ("R3 N29+LAI",  θN29, base_lai),
        ("R4 B+LAI",    θB,   base_lai),
        (@sprintf("R5 N29+pipe%+.0fK", PIPE_TEST_K), θN29, base_pipe)]
println("Running $(length(runs)) full-season forward simulations ...")
T2s = Vector{Vector{Float64}}()
for (name, θ, b) in runs
    t = @elapsed (T2 = collect(forward_map(θ, QoI_V, b, days, tspan, teval; u0 = copy(u0))[2]))
    length(T2) == NCOL || error("$name: solver stopped early ($(length(T2)) of $NCOL points)")
    push!(T2s, T2); @printf("  %-16s %.1f s\n", name, t)
end

# ── grid-level masks ─────────────────────────────────────────────────────────
coldays = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
okT   = .!isnan.(Tobs_g); okT[1] = false
night = A(Ig_v).(predtimes) .< 1.0
srcg  = A(Float64.(src)).(predtimes)
Uscr  = A(U1v).(predtimes)
Toutg = A(Tout_v).(predtimes) .- 273.15
mb(T, m) = any(m) ? mean(T[m] .- Tobs_g[m]) : NaN

println("\n", "="^92)
println("  Overall bias / RMSE [K]: calibration days 0-41 vs holdout 42-114")
println("="^92)
mc = okT .& (coldays .< 42); mh = okT .& (coldays .>= 42)
@printf("  %-16s %8s %8s %8s %8s %10s %10s\n", "run", "bias cal", "RMSE cal", "bias hol", "RMSE hol", "night hol", "day hol")
for (k, (name, _, _)) in enumerate(runs)
    T = T2s[k]
    @printf("  %-16s %+8.2f %8.2f %+8.2f %8.2f %+10.2f %+10.2f\n", name, mb(T, mc), sqrt(mean((T[mc] .- Tobs_g[mc]).^2)),
            mb(T, mh), sqrt(mean((T[mh] .- Tobs_g[mh]).^2)), mb(T, mh .& night), mb(T, mh .& .!night))
end

println("\n  Weekly NIGHT bias [K] (Iglob < 1 W/m2) and drivers")
@printf("  %-7s %6s %7s %7s %7s |", "days", "Tout", "heatEst", "pipeMs", "screen")
for (name, _, _) in runs; @printf(" %9s", split(name)[1]); end
@printf(" | %8s\n", "R2 day")
for w0 in 0:7:N_DAYS-1
    wm = okT .& (coldays .>= w0) .& (coldays .< min(w0 + 7, N_DAYS)); wn = wm .& night
    any(wn) || continue
    @printf("  %3d-%-3d %6.1f %6.0f%% %6.0f%% %7.2f |", w0, min(w0 + 6, N_DAYS - 1), mean(Toutg[wn]),
            100mean(srcg[wn] .== 2), 100mean(srcg[wn] .== 1), mean(Uscr[wn]))
    for T in T2s; @printf(" %+9.2f", mb(T, wn)); end
    @printf(" | %+8.2f\n", mb(T2s[2], wm .& .!night))
end

println("\n  Holdout NIGHT bias by Tpipe source [K]  (n = grid points)")
@printf("  %-16s %18s %18s %18s\n", "run", "measured pipe", "heating estimate", "cold fallback")
for (k, (name, _, _)) in enumerate(runs)
    @printf("  %-16s", name)
    for s in 1:3
        m = mh .& night .& (srcg .== s); @printf(" %+9.2f (n=%5d)", mb(T2s[k], m), count(m))
    end
    println()
end

# daily night-bias correlations (R2)
dn  = [mb(T2s[2], okT .& night .& (coldays .== d)) for d in 0:N_DAYS-1]
dTo = [let m = okT .& night .& (coldays .== d); any(m) ? mean(Toutg[m]) : NaN end for d in 0:N_DAYS-1]
dHe = [let m = okT .& night .& (coldays .== d); any(m) ? mean(srcg[m] .== 2) : NaN end for d in 0:N_DAYS-1]
dSc = [let m = okT .& night .& (coldays .== d); any(m) ? mean(Uscr[m]) : NaN end for d in 0:N_DAYS-1]
ok = .!isnan.(dn) .& .!isnan.(dTo)
@printf("\n  Daily night bias (R2) correlation:  vs Tout %+.2f   vs heat-estimate share %+.2f   vs screen %+.2f\n",
        cor(dn[ok], dTo[ok]), cor(dn[ok], dHe[ok]), cor(dn[ok], dSc[ok]))

# ── figures ──────────────────────────────────────────────────────────────────
dd = 0:N_DAYS-1
p1 = plot(xlabel = "day", ylabel = "night bias [K]", title = "Daily NIGHT bias (model − measured)", legend = :bottomleft, size = (1100, 380))
for (k, (name, _, _)) in enumerate(runs)
    plot!(p1, dd, [mb(T2s[k], okT .& night .& (coldays .== d)) for d in dd], lw = 2, label = name)
end
hline!(p1, [0], lc = :black, label = ""); vline!(p1, [41.5], lc = :blue, ls = :dash, label = "end calibration")
p2 = plot(dd, 100 .* dHe, lc = :firebrick, lw = 2, label = "heating-estimate share of night [%]", xlabel = "day",
          ylabel = "%", legend = :topleft, size = (1100, 260))
plot!(p2, dd, 100 .* dSc, lc = :teal, lw = 2, label = "night screen U1 [%]")
p3 = plot(dd, 2.0 .* ones(length(dd)), lc = :gray, lw = 2, ls = :dash, label = "LAI_FIX (Vanthoor)", xlabel = "day",
          ylabel = "LAI", legend = :bottomright, size = (1100, 220))
plot!(p3, dd, lai_of_day.(dd), lc = :green, lw = 2, label = "derived observed LAI (R3/R4)")
savefig(plot(p1, p2, p3, layout = grid(3, 1, heights = [0.5, 0.3, 0.2]), size = (1100, 900)), DG_FIG * "_daily.png")
p4 = scatter(dTo[ok], dn[ok], zcolor = dHe[ok], m = (:circle, 5), c = :viridis, xlabel = "night mean Tout [°C]",
             ylabel = "night bias R2 [K]", colorbar_title = "heat-est share", title = "Night bias vs outdoor T (colour = heating-estimate share)",
             label = "", size = (700, 500))
savefig(p4, DG_FIG * "_scatter.png")
println("\nFigures: ", DG_FIG, "_daily.png, ", DG_FIG, "_scatter.png")
