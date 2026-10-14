# =============================================================================
# diagnose_vent_mapping.jl -- why the n54 physics needs its own vent controller (Note 57; forward only)
# =============================================================================
# Run from the GreenhouseSim folder (~3-4 min, no MCMC):
#     PRIORS=n54 julia --project=. test/diagnose_vent_mapping.jl
#
# At the n54 posterior, closed loop (CL) fits by day while open loop with the logged vents (OL) is ~2 K too warm
# and ~4 hPa too humid, at nearly the same mean opening (Note 56). Two parts:
#  A. Opening by hour of day: CL effective opening (max(u7, logged VentMin)) vs the logged (lee + wind)/200,
#     and lee / wind sides separately; for days 0-48 and 49-114.
#  B. Open-loop variants of the vent INPUT (physics unchanged):
#       OL mean    (VentLee + Ventwind)/200          [the current mapping]
#       OL max     max(VentLee, Ventwind)/100        [both sides as wide as the wider one]
#       OL lee     VentLee/100 ;  OL wind  Ventwind/100
#       OL x1.5    1.5 x mean (clamped)              [how much more opening the physics needs]
#       OL CLrep   the CL effective opening replayed  [sanity: should reproduce CL]
#     If one mapping fixes the day (T and VP) the issue is how the log is read; if only a scale factor does, the
#     vent flow per opening is too low; if only CLrep does, it is the TIMING of the opening.
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n54")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2")
    isfile(f) || error("missing $f")
    push!(parts, JLD2.load(f)["chain"][30001:end, 1:NM])
end
pm = vec(mean(reduce(vcat, parts), dims = 1)); pP = update_params(params_cl, CQ_MODEL, pm)
pOL = update_params(pP, ["vent_ctrl"], [0.0])

# raw lee / wind logs on the same rows as ts
Leev = Float64[]; Windv = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i])
    push!(Leev, isnan(vl) ? 0.0 : clamp(vl / 100, 0.0, 1.0)); push!(Windv, isnan(vw) ? 0.0 : clamp(vw / 100, 0.0, 1.0))
end
length(Leev) == length(ts) || error("row filter mismatch")

ctrl_with(U8vec) = ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, zc, Z(U8vec), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
                                  Z(VentSPv), Z(HumSPv), Z(Blackv), Z(VentMinv))
function simdays_c(p, ds, u0x, rtol, ctrl)
    us = Matrix{Float64}(undef, length(u0x), 1 + length(ds) * NPD); us[:, 1] = u0x
    u = copy(u0x); col = 2
    for d in ds
        ctx = GreenhouseSim.SimContext(p, weather_c, ctrl, t0 + 86400.0 * d, 86400.0 * d)
        sol = try
            GreenhouseSim.solve(GreenhouseSim.ODEProblem(GreenhouseSim.rhs!, u, tspan, ctx), GreenhouseSim.Tsit5();
                                saveat = teval, reltol = rtol)
        catch
            return nothing
        end
        length(sol.t) == NPD || return nothing
        S = Array(sol); us[:, col:col+NPD-1] = S; col += NPD; u = S[:, end]
    end
    return us
end

okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g); okU = .!isnan.(U8obs); okU[1] = false
inset(ds) = in.(cd_, Ref(Set(ds)))
VP(T, RH) = RH ./ 100 .* Pws.(T) ./ 100; VPo = VP(Tobs_g, RHobs_g)
mb(v, o, m) = any(m) ? mean(v[m] .- o[m]) : NaN
mm(v, m) = any(m) ? mean(v[m]) : NaN
HB = Dict{Tuple{Int,Int},Vector{Int}}()
for c in 2:NCOL; push!(get!(HB, (cd_[c], floor(Int, ((c - 2) % NPD) * 300 / 3600)), Int[]), c); end
function hourly_rmse(v, o, mask)
    e = Float64[]
    for cols in values(HB); k = [c for c in cols if mask[c]]; length(k) >= 6 && push!(e, mean(v[k]) - mean(o[k])); end
    isempty(e) ? NaN : sqrt(mean(abs2, e))
end
PER = [(0, 8), (9, 48), (49, 100), (101, 109), (110, 114)]
hrc = [floor(Int, mod(((c - 2) % NPD) * 300 / 3600, 24.0)) for c in 1:NCOL]; hrc[1] = 0
HB_ = ((6, 9), (9, 12), (12, 15), (15, 18), (18, 21))

# ---- CL run and its effective opening -------------------------------------------------------------------------
usCL = simdays(pP, days, u0_cl, RELTOL); usCL === nothing && error("CL solve failed")
VminC = A(VentMinv).(predtimes)
U8cl = max.(clamp.(usCL[7, :], 0.0, 1.0), VminC)
LeeC = A(Leev).(predtimes); WindC = A(Windv).(predtimes)

println("\n", "="^120)
println("  A. Vent opening by hour of day (hourly means of 5-min values; fraction 0-1)")
for (lab, dsel) in (("days 0-48", collect(0:48)), ("days 49-114", collect(49:114)))
    md = inset(dsel) .& okU
    @printf("\n  %-11s | %s\n", lab, join([@sprintf("%5d", h) for h in 0:2:22], " "))
    for (nm_, v) in (("CL model", U8cl), ("log mean", U8obs), ("log lee", LeeC), ("log wind", WindC), ("VentMin", VminC))
        @printf("  %-11s | %s\n", nm_, join([@sprintf("%5.2f", mm(v, md .& (hrc .>= h) .& (hrc .< h + 2))) for h in 0:2:22], " "))
    end
    @printf("  %-11s | %s\n", "CL - log", join([@sprintf("%+5.2f", mm(U8cl .- U8obs, md .& (hrc .>= h) .& (hrc .< h + 2))) for h in 0:2:22], " "))
end
@printf("\n  correlation of 5-min openings, CL vs log: day %.2f  night %.2f;   lee vs wind (log): %.2f\n",
        cor(U8cl[okU .& dayt], U8obs[okU .& dayt]), cor(U8cl[okU .& nite], U8obs[okU .& nite]), cor(LeeC[okU], WindC[okU]))

# ---- B. OL variants ---------------------------------------------------------------------------------------------
U8clts = stepfun(predtimes, U8cl).(ts)
variants = [("CL (n54 as fitted)", pP, ctrl_x),
            ("OL mean (current)", pOL, ctrl_x),
            ("OL max(lee,wind)", pOL, ctrl_with(max.(Leev, Windv))),
            ("OL lee only", pOL, ctrl_with(Leev)),
            ("OL wind only", pOL, ctrl_with(Windv)),
            ("OL 1.5 x mean", pOL, ctrl_with(clamp.(1.5 .* U8z, 0.0, 1.0))),
            ("OL 2 x mean", pOL, ctrl_with(clamp.(2.0 .* U8z, 0.0, 1.0))),
            ("OL CL replayed", pOL, ctrl_with(U8clts))]
# Note 57: concave opening -> flow (more flow per % at small openings), logged mean mapping, physics otherwise n54
for ve in (0.8, 0.65, 0.5)
    push!(variants, ("OL mean, vent_exp $(ve)", update_params(pOL, ["vent_exp"], [ve]), ctrl_x))
end
push!(variants, ("OL lee, vent_exp 0.65", update_params(pOL, ["vent_exp"], [0.65]), ctrl_with(Leev)))
ONLY_EXP = get(ENV, "ONLY_EXP", "0") == "1"      # ONLY_EXP=1: skip part A and the first B rows (rerun after Note 57)
ONLY_EXP && (variants = [v for v in variants if occursin("vent_exp", v[1]) || startswith(v[1], "OL mean (")])
println("\n", "="^170)
@printf("  B. %-20s | %6s | %s | %s | %6s %6s %6s | %s\n", "run", "hRMS", join([@sprintf("%6s", "n$(a)") for (a, b) in PER], " "),
        join([@sprintf("%6s", "d$(a)-$(b)") for (a, b) in HB_], " "), "T day", "VPday", "VPnit", "mean opening day / night")
println("  ", "-"^168)
for (nm_, p, ctrl) in variants
    us = simdays_c(p, days, u0_cl, RELTOL, ctrl)
    us === nothing && (println("  $nm_: solver failed"); continue)
    T = us[2, :]; RH = rhf.(us[2, :], us[4, :])
    Uin = startswith(nm_, "CL") ? U8cl : A(ctrl === ctrl_x ? U8z : ctrl.U8.(trel)).(predtimes)
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    @printf("  %-23s | %6.3f | %s | %s | %+6.2f %+6.2f %+6.2f | %.3f / %.3f\n", nm_, hourly_rmse(T, Tobs_g, okT),
            join([@sprintf("%+6.2f", v) for v in nb], " "), join([@sprintf("%+6.2f", v) for v in hb], " "),
            mb(T, Tobs_g, okT .& dayt), mb(VP(T, RH), VPo, okV .& dayt), mb(VP(T, RH), VPo, okV .& nite),
            mm(Uin, okU .& dayt), mm(Uin, okU .& nite))
end
println("""
  Bias = model - measured (K; VP in hPa). All runs: n54 posterior mean, full season from day 0, AiCU inputs.
  Read: which OL input reproduces the CL day (T and VP)? A mapping (max / lee / wind) -> the log is read wrongly;
  only a scale factor -> vent flow per opening too low; only CL replayed -> timing of the opening matters.""")
