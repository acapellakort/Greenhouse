# =============================================================================
# diagnose_dusk_events.jl -- why does the model crash when the grower vents at dusk / night? (Note 39)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~1 min, no MCMC):
#     julia --project=. test/diagnose_dusk_events.jl
#
# Event = onset of a logged vent opening at dusk/night: U8_log >= 0.15 while the previous 30 min averaged
# < 0.10, Iglobal < 50 W/m2, at least 2 h after the previous event. Whole season (days 1-114).
# Response = mean T over +45..+75 min minus mean T over -15..0 min (measured and model).
# Models, all at the Note 38 posterior mean (PRIORS=n37):
#   M_cl     closed-loop vents (Note 38 as calibrated)                    -- reference
#   M_log    logged vents replayed ALL day (vc_night_log=1, vc_night_I huge): the model sees the same openings
#   M_log_s  as M_log with vent_scr_perm = 0.2 (vents strongly blocked by a closed screen)
# Screen hypothesis: if the measured response is weak when the screen is closed (U1 >= 0.5) and M_log
# over-reacts mainly then, while M_log_s fixes those events, a screen-dependent night vent efficiency is the
# missing structure. Heating hypothesis: if M_log falls >1 K below the heating setpoint after the event, the
# PI heating (capped at hc_Tmax) cannot compensate.
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n37"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))
for k in (:vc_night_log, :vc_min_night, :vc_night_I)
    hasproperty(params, k) || error("constants_climate.json has no $k entry (Note 38/39 flags missing)")
end

parts = Vector{Vector{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n37_seed$(s).jld2")
    isfile(f) || error("missing $f (Note 38 chain)")
    d = JLD2.load(f); push!(parts, vec(mean(d["chain"][d["burn_in"]+1:end, 1:NM], dims = 1)))
end
pm = (parts[1] .+ parts[2]) ./ 2
function run_p(names, vals)
    p = update_params(params_cl, CQ_MODEL, pm)
    isempty(names) || (p = update_params(p, names, vals))
    us = simdays(p, days, u0_8, RELTOL); us === nothing && error("solver failed for $(names)")
    us[2, :]
end
models = [("M_cl",    String[], Float64[]),
          ("M_log",   ["vc_night_log", "vc_night_I"], [1.0, 1e9]),
          ("M_log_s", ["vc_night_log", "vc_night_I", "vent_scr_perm"], [1.0, 1e9, 0.2])]
println("\nRunning $(length(models)) full-season simulations ...")
Tm = Dict{String,Vector{Float64}}()
for (nm_, nn, vv) in models; t = @elapsed (Tm[nm_] = run_p(nn, vv)); @printf("  %-8s %.1f s\n", nm_, t); end

# ── time helpers on the 5-min model grid ─────────────────────────────────────
function colat(t)
    x = (t - t0) / 86400; d = floor(Int, x); j = round(Int, (x - d) * 86400 / 300)
    j >= NPD - 1 && (d += 1; j = 0)
    (d < 0 || d > last(days)) ? 0 : col_of(d, j)
end
wmean(v, te, k1, k2) = (cs = [colat(te + 300k) for k in k1:k2]; cs = [c for c in cs if c > 0 && !isnan(v[c])];
                        length(cs) >= (k2 - k1 + 1) ÷ 2 ? mean(v[cs]) : NaN)
AU1 = A(U1v); AW = A(Wind_v); ATo = A(Tout_v); AU11 = A(U11v); APp = A(Tpipe_meas); AIg = A(Ig_v)

# ── event detection (on the logged vent, 5-min grid) ─────────────────────────
ev = Float64[]; last_te = -Inf
for d in 1:last(days), j in 0:NPD-2
    te = t0 + 86400.0 * d + 300.0 * j
    u = AU(te); (isnan(u) || u < 0.15) && continue
    AIg(te) < 50 || continue
    pre = [AU(te - 300k) for k in 1:6]; any(isnan, pre) && continue
    mean(pre) < 0.10 || continue
    te - last_te >= 7200 || continue
    push!(ev, te); global last_te = te
end
@printf("\n%d dusk/night vent-opening events found (U8 >= 0.15 after < 0.10, Iglobal < 50 W/m2)\n", length(ev))

rows = NamedTuple[]
for te in ev
    Tob = Tobs_g
    dTo = wmean(Tob, te, 9, 15) - wmean(Tob, te, -3, 0)
    isnan(dTo) && continue
    sv = filter(!isnan, [AU1(te + 300k) for k in 0:6]); uv = filter(!isnan, [AU(te + 300k) for k in 0:12])
    (isempty(sv) || isempty(uv)) && continue
    r = (te = te, day = floor(Int, (te - t0) / 86400), hour = ((te - t0) % 86400) / 3600,
         scr = mean(sv), umax = maximum(uv), umean = mean(uv),
         wind = AW(te), dTout = AT(te) - ATo(te), pipe = APp(te) - 273.15, sp = AU11(te),
         dTo = dTo, dTcl = wmean(Tm["M_cl"], te, 9, 15) - wmean(Tm["M_cl"], te, -3, 0),
         dTlog = wmean(Tm["M_log"], te, 9, 15) - wmean(Tm["M_log"], te, -3, 0),
         dTlogs = wmean(Tm["M_log_s"], te, 9, 15) - wmean(Tm["M_log_s"], te, -3, 0),
         Tlog1h = wmean(Tm["M_log"], te, 9, 15))
    push!(rows, r)
end

# ── summary by screen state ──────────────────────────────────────────────────
med(x) = (y = filter(!isnan, collect(x)); isempty(y) ? NaN : median(y))
mabs(x) = (y = filter(!isnan, collect(x)); isempty(y) ? NaN : mean(abs.(y)))
function slope(x, y)   # least-squares slope of y on x through the origin
    k = [i for i in eachindex(x) if !isnan(x[i]) && !isnan(y[i])]
    isempty(k) ? NaN : sum(x[k] .* y[k]) / sum(x[k] .^ 2)
end
groups = [("screen closed (U1 >= 0.5)", r -> r.scr >= 0.5), ("screen open (U1 < 0.5)", r -> r.scr < 0.5), ("all", r -> true)]
println("\n", "="^128)
println("  Response to a dusk/night vent opening: ΔT = T(+45..75 min) - T(-15..0 min), K. 'sens' = median ΔT / (mean U8 x (Tair-Tout)).")
println("  'crash' = model drop > 2 K while measured drop < 0.5 K. 'heat-sat' = M_log ends > 1 K below the heating setpoint.")
@printf("  %-26s %4s | %7s %7s %7s %7s | %7s %7s %7s | %6s %6s | %6s %6s | %8s\n", "group", "n", "ΔT meas", "ΔT cl",
        "ΔT log", "ΔT logs", "|err|cl", "|err|lg", "|err|ls", "sensM", "sensL", "crashL", "crashS", "heat-sat")
println("  ", "-"^126)
for (gname, f) in groups
    g = filter(f, rows); n = length(g); n == 0 && continue
    sens(v) = med([v[i] / (g[i].umean * max(g[i].dTout, 1.0)) for i in eachindex(g)])
    crash(v) = count(i -> v[i] < -2 && g[i].dTo > -0.5, eachindex(g))
    @printf("  %-26s %4d | %+7.2f %+7.2f %+7.2f %+7.2f | %7.2f %7.2f %7.2f | %+6.3f %+6.3f | %6d %6d | %8d\n", gname, n,
            med(r.dTo for r in g), med(r.dTcl for r in g), med(r.dTlog for r in g), med(r.dTlogs for r in g),
            mabs(r.dTcl - r.dTo for r in g), mabs(r.dTlog - r.dTo for r in g), mabs(r.dTlogs - r.dTo for r in g),
            sens([r.dTo for r in g]), sens([r.dTlog for r in g]),
            crash([r.dTlog for r in g]), crash([r.dTlogs for r in g]),
            count(r -> r.Tlog1h < r.sp - 1.0, g))
end
gc = filter(r -> r.scr >= 0.5, rows); go = filter(r -> r.scr < 0.5, rows)
@printf("\n  Slope of model ΔT on measured ΔT (1 = right sensitivity): screen closed  M_log %.2f  M_log_s %.2f  |  screen open  M_log %.2f  M_log_s %.2f\n",
        slope([r.dTo for r in gc], [r.dTlog for r in gc]), slope([r.dTo for r in gc], [r.dTlogs for r in gc]),
        slope([r.dTo for r in go], [r.dTlog for r in go]), slope([r.dTo for r in go], [r.dTlogs for r in go]))
@printf("  Event conditions (median): screen-closed events  wind %.1f m/s, Tair-Tout %.1f K, pipe %.1f °C, peak U8 %.2f, hour %.1f\n",
        med(r.wind for r in gc), med(r.dTout for r in gc), med(r.pipe for r in gc), med(r.umax for r in gc), med(r.hour for r in gc))
@printf("                              screen-open events    wind %.1f m/s, Tair-Tout %.1f K, pipe %.1f °C, peak U8 %.2f, hour %.1f\n",
        med(r.wind for r in go), med(r.dTout for r in go), med(r.pipe for r in go), med(r.umax for r in go), med(r.hour for r in go))

# worst 12 M_log over-reactions
println("\n  12 largest M_log over-reactions (model ΔT - measured ΔT):")
@printf("  %4s %5s %5s %5s %5s %5s %6s | %7s %7s %7s %7s\n", "day", "hour", "scr", "U8max", "wind", "dT_o", "pipe", "meas", "M_cl", "M_log", "M_log_s")
for r in sort(rows, by = r -> r.dTlog - r.dTo)[1:min(12, length(rows))]
    @printf("  %4d %5.1f %5.2f %5.2f %5.1f %5.1f %6.1f | %+7.2f %+7.2f %+7.2f %+7.2f\n", r.day, r.hour, r.scr, r.umax, r.wind,
            r.dTout, r.pipe, r.dTo, r.dTcl, r.dTlog, r.dTlogs)
end
println("""
  How to read it:
   - |err| and crashes concentrated in 'screen closed', and M_log_s fixes them -> screen-dependent night vent
     efficiency (default-off flag next, then calibrate).
   - crashes with heat-sat > 0 and similar in both screen states -> heating capacity / pipe saturation.
   - M_log over-reacts in both screen states similarly -> general vent sensitivity at low Tair-Tout or low wind.""")

# ── figure: measured vs model response ───────────────────────────────────────
pl = plot(xlabel = "measured ΔT (K)", ylabel = "model ΔT (K)", title = "Dusk/night vent events: 1-h response",
          legend = :bottomright, size = (700, 600))
plot!(pl, [-6, 3], [-6, 3], lc = :gray, ls = :dash, label = "1:1")
scatter!(pl, [r.dTo for r in gc], [r.dTlog for r in gc], ms = 4, label = "M_log, screen closed")
scatter!(pl, [r.dTo for r in go], [r.dTlog for r in go], ms = 4, label = "M_log, screen open")
scatter!(pl, [r.dTo for r in gc], [r.dTlogs for r in gc], ms = 3, marker = :x, label = "M_log_s, screen closed")
fig = joinpath(FIGDIR, "AiCU_dusk_events.png"); savefig(pl, fig); println("Figure: ", fig)
