# =============================================================================
# inspect_spike_days.jl -- what is different on the spike days (24, 30, 41)?
# =============================================================================
# Run from the GreenhouseSim folder (~30 s, no MCMC, no src/ changes):
#     julia --project=. test/inspect_spike_days.jl
#     SHOW_DAYS=21,24,30,38,41 julia --project=. test/inspect_spike_days.jl
#
# Two forward runs over days 0-41 (Note 28/29 climate values, U8_min_leak = 0 and 0.05),
# then a per-DAYTIME table (Iglob > 50 W/m2) for every day, built from the RAW logger
# columns, and 5-min panels for the selected days. Questions it answers:
#   1. Are the logged vents really ~closed while it is hot inside (vs the VIP vent setpoint)?
#   2. Are vent/temperature readings missing (NaN -> mapped to 0 = closed in the model)?
#   3. Is the BLACKOUT/SHADE screen (BlackScr) closed during the day? The model ignores
#      BlackScr entirely (only EnScr -> U1), so a shaded day gets full sun in the model.
# =============================================================================

ENV["OPTB_MODE"] = "setup"
ENV["N_DAYS"]    = "42"
include(joinpath(@__DIR__, "run_joint_inference_optionB.jl"))

const SHOW = parse.(Int, split(get(ENV, "SHOW_DAYS", "21,24,30,38,41"), ","))
const SP_FIG = joinpath(FIGDIR, "AiCU_spike_days")

# ── raw logger rows in the window (same row filter as the joint script) ─────
raw = Dict(k => Float64[] for k in ("Tair", "VentLee", "Ventwind", "EnScr", "BlackScr", "Iglob", "Tout",
                                     "VentSP", "HeatSP", "PipeGrow", "RHair"))
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    push!(raw["Tair"], _num(gh.Tair[i]));       push!(raw["VentLee"], _num(gh.VentLee[i]))
    push!(raw["Ventwind"], _num(gh.Ventwind[i])); push!(raw["EnScr"], _num(gh.EnScr[i]))
    push!(raw["BlackScr"], _num(gh.BlackScr[i])); push!(raw["Iglob"], _num(mt.Iglob[i]))
    push!(raw["Tout"], _num(mt.Tout[i]));        push!(raw["PipeGrow"], _num(gh.PipeGrow[i]))
    push!(raw["RHair"], _num(gh.RHair[i]))
    vs = filter(!isnan, vcat([_num(vip[i, Symbol("VentLeew_Vip_$z")]) for z in 1:6],
                             [_num(vip[i, Symbol("VentWind_Vip_$z")]) for z in 1:6]))
    push!(raw["VentSP"], isempty(vs) ? NaN : mean(vs))
    hs = filter(!isnan, [_num(vip[i, Symbol("HeatTemp_Vip_$z")]) for z in (1, 2, 3, 5, 6)])
    push!(raw["HeatSP"], isempty(hs) ? NaN : mean(hs))
end
length(raw["Tair"]) == length(ts) || error("row filter mismatch")
rday  = floor.(Int, (ts .- t0) ./ 86400)
rhour = mod.((ts .- t0) ./ 3600, 24)

# ── model runs ───────────────────────────────────────────────────────────────
θN = [4.46e-5, 44.0, 24.5, 0.05]
T_l05 = collect(forward_map(θN, QoI_V, base, days, tspan, teval; u0 = copy(u0))[2]) .- 273.15
T_l00 = collect(forward_map([θN[1:3]; 0.0], QoI_V, base, days, tspan, teval; u0 = copy(u0))[2]) .- 273.15
# model value at each logger row (nearest grid point of that day)
gcol(i) = col_of(rday[i], clamp(round(Int, (ts[i] - t0 - 86400.0 * rday[i]) / 300.0), 0, NPD - 1))
Mrow05 = [T_l05[gcol(i)] for i in eachindex(ts)]
Mrow00 = [T_l00[gcol(i)] for i in eachindex(ts)]

nanf(v) = isempty(v) ? NaN : mean(isnan.(v))
mx(v) = (w = filter(!isnan, v); isempty(w) ? NaN : maximum(w))
mn(v) = (w = filter(!isnan, v); isempty(w) ? NaN : mean(w))

println("\n", "="^128)
println("  DAYTIME (Iglob > 50 W/m²) per day, days 0-41.  U8 = mean(VentLee, Ventwind)/100 as the model sees it (NaN -> 0)")
println("="^128)
@printf("  %3s %6s %6s %6s %6s | %6s %6s %6s %6s %6s | %6s %6s | %6s %6s %6s %6s | %s\n",
        "day", "Iglob", "Tout", "Tmeas", "Tmeas", "M l=0", "M l=.05", "VentSP", "Tm-SP", "Tm>SP", "U8", "NaN", "EnScr", "Black", "RH", "pipe", "flag")
@printf("  %3s %6s %6s %6s %6s | %6s %6s %6s %6s %6s | %6s %6s | %6s %6s %6s %6s |\n",
        "", "max", "mean", "mean", "max", "max", "max", "mean", "max", "%time", "mean", "vent", "mean", "mean", "mean", "max")
for d in 0:N_DAYS-1
    m = (rday .== d) .& .!isnan.(raw["Iglob"]) .& (raw["Iglob"] .> 50)
    any(m) || continue
    ventraw = (raw["VentLee"][m] .+ raw["Ventwind"][m]) ./ 2
    u8 = [isnan(v) ? 0.0 : v / 100 for v in ventraw]
    above = raw["Tair"][m] .- raw["VentSP"][m]
    flag = (mx(Mrow00[m]) - mx(raw["Tair"][m]) > 8) ? "SPIKE (model max > meas max + 8 K)" : ""
    @printf("  %3d %6.0f %6.1f %6.1f %6.1f | %6.1f %6.1f %6.1f %+6.1f %5.0f%% | %6.3f %5.0f%% | %6.0f %6.0f %6.0f %6.1f | %s\n",
            d, mx(raw["Iglob"][m]), mn(raw["Tout"][m]), mn(raw["Tair"][m]), mx(raw["Tair"][m]),
            mx(Mrow00[m]), mx(Mrow05[m]), mn(raw["VentSP"][m]), mx(above), 100 * mn(Float64.(above .> 0)),
            mean(u8), 100 * nanf(ventraw), mn(raw["EnScr"][m]), mn(raw["BlackScr"][m]), mn(raw["RHair"][m]),
            mx(raw["PipeGrow"][m]), flag)
end
println("""
  Read: on a real spike day with logged vents ~0, the measured air should sit far above the VIP vent
  setpoint (Tm-SP >> 0) — if instead Tmeas stays near the setpoint while U8 ~ 0, the greenhouse was being
  cooled by something the model does not see (unlogged venting, shade/blackout screen, misting).
  Black = blackout/shade screen %: NOT an input to the model.""")

# ── 5-min panels for selected days ───────────────────────────────────────────
cols = Any[]
for d in SHOW
    m = rday .== d; h = rhour[m]
    pT = plot(h, raw["Tair"][m], lc = :black, lw = 2, label = "Tair meas", title = "day $d", titlefontsize = 10,
              ylabel = d == SHOW[1] ? "°C" : "", legend = d == SHOW[1] ? :topleft : false, legendfontsize = 6)
    plot!(pT, h, Mrow00[m], lc = :darkorange, lw = 1.5, label = "model leak=0")
    plot!(pT, h, Mrow05[m], lc = :purple, lw = 1.5, ls = :dash, label = "model leak=0.05")
    plot!(pT, h, raw["VentSP"][m], lc = :red, ls = :dot, label = "VIP vent SP")
    plot!(pT, h, raw["Tout"][m], lc = :steelblue, lw = 1, label = "Tout")
    pV = plot(h, raw["VentLee"][m], lc = :teal, label = "VentLee %", ylims = (-5, 105),
              ylabel = d == SHOW[1] ? "%" : "", legend = d == SHOW[1] ? :topleft : false, legendfontsize = 6)
    plot!(pV, h, raw["Ventwind"][m], lc = :green, ls = :dash, label = "Ventwind %")
    pS = plot(h, raw["EnScr"][m], lc = :brown, label = "EnScr % (model U1)", ylims = (-5, 105),
              legend = d == SHOW[1] ? :topleft : false, legendfontsize = 6)
    plot!(pS, h, raw["BlackScr"][m], lc = :black, ls = :dash, label = "BlackScr % (NOT in model)")
    pI = plot(h, raw["Iglob"][m], lc = :goldenrod, label = "Iglob W/m²", xlabel = "hour",
              legend = d == SHOW[1] ? :topleft : false, legendfontsize = 6)
    append!(cols, [pT, pV, pS, pI])
end
nd = length(SHOW)
order = [cols[4 * (c - 1) + r] for r in 1:4 for c in 1:nd]     # row-major: T row, vent row, screen row, light row
savefig(plot(order..., layout = (4, nd), size = (330 * nd, 900)), SP_FIG * ".png")
println("\nFigure: ", SP_FIG, ".png")
