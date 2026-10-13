# =============================================================================
# diagnose_night_heat_balance.jl -- which heat-balance term causes the night cold bias?
# =============================================================================
# Run from the GreenhouseSim folder (~1 min, no MCMC, no src/ changes):
#     julia --project=. test/diagnose_night_heat_balance.jl
#
# All runs start from the Note 28/29 climate values (N29) and change ONE existing
# parameter through update_params (no engine edits):
#   S0 baseline           N29
#   S1 eps_screen = 0     removes r12, the screen->sky longwave term that is ADDED to the
#                         air's losses when the screen closes (climate.jl r12; it was added
#                         for the venting calibration, where U1 = 0, and never tested at night)
#   S2 n_pipes 2 -> 3     +50% pipe->air convection (h4)
#   S3 k_ground = 50      the joint-posterior value: floor anchored to a FIXED 18 C deep soil
#   S4 S1 + pipe +10 K    S1 plus the +10 K heating-estimate test from diagnose_heating_bias.jl
# Then splits the night bias by screen state, and regresses the S0 night residual on
# (T_in - T_out), screen, and (T_pipe - T_in) on measured-pipe nights.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = get(ENV, "N_FULL", "115")
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const NH_FIG = joinpath(FIGDIR, "AiCU_night_heat_balance")
@printf("\nEngine values: eps_screen=%.2f epsil6=%.2f tau3=%.2f n_pipes=%.0f phi1=%.3f gamma1=%.2f\n",
        params.eps_screen, params.epsil6, params.tau3, params.n_pipes, params.phi1, params.gamma1)

# Tpipe source per row (same rules as the joint script), for the +10 K test and the regression
srcv = Int[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    To = _num(mt.Tout[i]); To_val = isnan(To) ? 15.0 : To
    pg = _num(gh.PipeGrow[i]); pl = _num(gh.PipeLow[i])
    hv = filter(!isnan, [_num(vip.HeatTemp_Vip_1[i]), _num(vip.HeatTemp_Vip_2[i]), _num(vip.HeatTemp_Vip_3[i]),
                         _num(vip.HeatTemp_Vip_5[i]), _num(vip.HeatTemp_Vip_6[i])])
    Th = isempty(hv) ? 20.0 : mean(hv)
    valid = (!isnan(pg) && pg > PIPE_MIN_VALID) || (!isnan(pl) && pl > PIPE_MIN_VALID)
    push!(srcv, valid ? 1 : (To_val < Th - 2.0 && To_val < PIPE_WARM_TOUT) ? 2 : 3)
end
length(srcv) == length(ts) || error("row filter mismatch")
Tpipe_alt = [srcv[i] == 2 ? Tpipev[i] + 10.0 : Tpipev[i] for i in eachindex(Tpipev)]
base_pipe = SimBase(params, weather,
                    ControlInputs(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8v), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_alt)), t0)

nV = QoI_V; θN = [4.46e-5, 44.0, 24.5, 0.05]
const RUNSET = get(ENV, "RUNSET", "A")
runs = RUNSET == "A" ? [
        ("S0 N29",            nV,                          θN,                               base),
        ("S1 eps_screen=0",   vcat(nV, "eps_screen"),      vcat(θN, 0.0),                    base),
        ("S2 n_pipes=3",      vcat(nV, "n_pipes"),         vcat(θN, 3.0),                    base),
        ("S3 k_ground=50",    nV,                          [θN[1], θN[2], 50.0, θN[4]],      base),
        ("S4 S1+pipe+10K",    vcat(nV, "eps_screen"),      vcat(θN, 0.0),                    base_pipe)] :
    [   # RUNSET B: which LOSS term scales with (T_in - T_out)?  one change at a time from N29
        ("S0 N29",            nV,                          θN,                               base),
        ("L1 U8_min_leak=0",  nV,                          [θN[1], θN[2], θN[3], 0.0],       base),
        ("L2 epsil6=0.22",    vcat(nV, "epsil6"),          vcat(θN, 0.22),                   base),
        ("L3 epsil6=0",       vcat(nV, "epsil6"),          vcat(θN, 0.0),                    base),
        ("L4 h_cover_out=15", vcat(nV, "h_cover_out"),     vcat(θN, 15.0),                   base),
        ("L5 nu4=1e-5",       nV,                          [1.0e-5, θN[2], θN[3], θN[4]],    base)]
T2s = Vector{Vector{Float64}}()
println("Running $(length(runs)) full-season forward simulations ...")
for (name, nm, x, b) in runs
    T2 = collect(forward_map(x, nm, b, days, tspan, teval; u0 = copy(u0))[2])
    length(T2) == NCOL || error("$name: solver stopped early")
    push!(T2s, T2)
end

coldays = vcat(-1, [(c - 2) ÷ NPD for c in 2:NCOL])
okT   = .!isnan.(Tobs_g); okT[1] = false
night = A(Ig_v).(predtimes) .< 1.0
scr   = A(U1v).(predtimes)
srcg  = A(Float64.(srcv)).(predtimes)
Toutg = A(Tout_v).(predtimes)
Tpg   = A(Tpipev).(predtimes)
mb(T, m) = any(m) ? mean(T[m] .- Tobs_g[m]) : NaN
mc = okT .& (coldays .< 42); mh = okT .& (coldays .>= 42)
sc_on = scr .>= 0.5

println("\n", "="^104)
@printf("  %-17s %8s %8s %8s %8s | %10s %10s | %11s %11s\n", "run", "bias cal", "RMSE cal", "bias hol", "RMSE hol",
        "night all", "day all", "night scr<", "night scr>")
@printf("  %-17s %8s %8s %8s %8s | %10s %10s | %11s %11s\n", "", "", "", "", "", "(0-114)", "(0-114)", "0.5", "=0.5")
println("  ", "-"^102)
for (k, (name, _, _, _)) in enumerate(runs)
    T = T2s[k]
    @printf("  %-17s %+8.2f %8.2f %+8.2f %8.2f | %+10.2f %+10.2f | %+11.2f %+11.2f\n", name,
            mb(T, mc), sqrt(mean((T[mc] .- Tobs_g[mc]).^2)), mb(T, mh), sqrt(mean((T[mh] .- Tobs_g[mh]).^2)),
            mb(T, okT .& night), mb(T, okT .& .!night), mb(T, okT .& night .& .!sc_on), mb(T, okT .& night .& sc_on))
end
@printf("  night points with screen >= 0.5: %d of %d\n", count(okT .& night .& sc_on), count(okT .& night))

println("\n  Weekly NIGHT bias [K]")
@printf("  %-7s %6s %6s |", "days", "Tout", "screen"); for r in runs; @printf(" %9s", split(r[1])[1]); end; println()
for w0 in 0:7:N_DAYS-1
    wn = okT .& night .& (coldays .>= w0) .& (coldays .< min(w0 + 7, N_DAYS)); any(wn) || continue
    @printf("  %3d-%-3d %6.1f %6.2f |", w0, min(w0 + 6, N_DAYS - 1), mean(Toutg[wn]) - 273.15, mean(scr[wn]))
    for T in T2s; @printf(" %+9.2f", mb(T, wn)); end; println()
end

# regression of each run's night residual on measured-pipe nights (pipe temperature known there).
# The right physical fix should FLATTEN the (T_in - T_out) slope, not just shift the intercept.
m = okT .& night .& (srcg .== 1)
X = hcat(ones(count(m)), Tobs_g[m] .- Toutg[m], scr[m], Tpg[m] .- Tobs_g[m])
println("\n  Night residual regression, measured-pipe nights (n=$(count(m))):  resid = a + b*(T_in-T_out) + c*screen + d*(T_pipe-T_in)")
@printf("  %-17s %8s %8s %8s %8s %6s\n", "run", "a [K]", "b", "c [K]", "d", "R²")
for (k, (name, _, _, _)) in enumerate(runs)
    y = T2s[k][m] .- Tobs_g[m]; β = X \ y; r2 = 1 - sum(abs2, y .- X * β) / sum(abs2, y .- mean(y))
    @printf("  %-17s %+8.2f %+8.3f %+8.2f %+8.3f %6.2f\n", name, β[1], β[2], β[3], β[4], r2)
end
println("  (b < 0: model loses too much heat to outside; c < 0: closing the screen makes the model too cold;")
println("   d > 0: model gains too little from the pipes)")

dd = 0:N_DAYS-1
p1 = plot(xlabel = "day", ylabel = "night bias [K]", title = "Daily NIGHT bias — single-term sensitivity (from N29)",
          legend = :bottomleft, size = (1100, 420))
for (k, (name, _, _, _)) in enumerate(runs)
    plot!(p1, dd, [mb(T2s[k], okT .& night .& (coldays .== d)) for d in dd], lw = 2, label = name)
end
hline!(p1, [0], lc = :black, label = ""); vline!(p1, [41.5], lc = :blue, ls = :dash, label = "")
savefig(p1, NH_FIG * "_daily_$(RUNSET).png")
println("\nFigure: ", NH_FIG, "_daily_$(RUNSET).png")
