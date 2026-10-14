# =============================================================================
# diagnose_vent_split.jl -- leeward and windward roof vents as two flows, open loop (Note 58; forward only)
# =============================================================================
# Run from the GreenhouseSim folder (~4-5 min, no MCMC):
#     PRIORS=n54 julia --project=. test/diagnose_vent_split.jl
#
# Physics = n54 posterior mean (vent_exp set to 1; curvature now vs_exp), vent_ctrl = 0, logged openings as inputs:
# U8 = VentLee/100, U7 = Ventwind/100, flag vent_split = 1.
#   flow = 1/2 [ lee^e sqrt(B + c_lee n3 v^2) + g_w wind^e sqrt(B + c_wind n3 v^2) ]
# Grid: windward discharge g_w (1, 2, 3), windward wind-pressure c_wind (1, 5), curvature e (1, 0.65).
# Reference rows: CL (n54 as fitted), OL with the current single-opening mapping, and split with all defaults
# (must equal the current mapping up to vent_exp 1.016 -> 1).
# Read: which combination brings the open-loop day (T, VP; early and late season) to the CL level without making
# the evenings (18-21 h) and nights cold? g_w / c_wind > 1 = windward side more effective per % opened;
# e < 1 = more flow per % at small openings.
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

Leev = Float64[]; Windv = Float64[]
for i in 1:n
    g = _num(gh.GHtime[i]); (isnan(g) || g < win_lo || g >= win_hi) && continue
    vl = _num(gh.VentLee[i]); vw = _num(gh.Ventwind[i])
    push!(Leev, isnan(vl) ? 0.0 : clamp(vl / 100, 0.0, 1.0)); push!(Windv, isnan(vw) ? 0.0 : clamp(vw / 100, 0.0, 1.0))
end
length(Leev) == length(ts) || error("row filter mismatch")
ctrl_split = ControlInputsX(Z(U1v), zc, zc, zc, zc, zc, Z(Windv), Z(Leev), zc, Z(U10v), Z(U11v), Z(U12v), Z(Tpipe_meas),
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

okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g)
inset(ds) = in.(cd_, Ref(Set(ds)))
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
hrc = [floor(Int, mod(((c - 2) % NPD) * 300 / 3600, 24.0)) for c in 1:NCOL]; hrc[1] = 0
HB_ = ((6, 9), (9, 12), (12, 15), (15, 18), (18, 21))
EARLY = inset(collect(0:48)); LATE = inset(collect(49:114))

psplit(gw, cw, e) = update_params(pOL, ["vent_split", "vent_exp", "vs_gwind", "vs_cwind", "vs_clee", "vs_exp"],
                                  [1.0, 1.0, gw, cw, 1.0, e])
runs = Any[("CL (n54 as fitted)", pP, ctrl_x), ("OL mean (current)", pOL, ctrl_x), ("split defaults", psplit(1.0, 1.0, 1.0), ctrl_split)]
for e in (1.0, 0.65), cw in (1.0, 5.0), gw in (1.0, 2.0, 3.0)
    (e == 1.0 && cw == 1.0 && gw == 1.0) && continue
    push!(runs, ("split g$(Int(gw)) c$(Int(cw)) e$(e)", psplit(gw, cw, e), ctrl_split))
end

println("\n", "="^190)
@printf("  %-22s | %6s | %s | %s | %6s %6s %6s | %6s %6s | %s\n", "run", "hRMS",
        join([@sprintf("%6s", "n$(a)") for (a, b) in PER], " "), join([@sprintf("%6s", "d$(a)-$(b)") for (a, b) in HB_], " "),
        "T day", "VPday", "VPnit", "day<49", "day49+", "")
println("  ", "-"^188)
for (nm_, p, ctrl) in runs
    us = simdays_c(p, days, u0_cl, RELTOL, ctrl)
    us === nothing && (println("  $nm_: solver failed"); continue)
    T = us[2, :]; RH = rhf.(us[2, :], us[4, :])
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    @printf("  %-22s | %6.3f | %s | %s | %+6.2f %+6.2f %+6.2f | %+6.2f %+6.2f |\n", nm_, hourly_rmse(T, Tobs_g, okT),
            join([@sprintf("%+6.2f", v) for v in nb], " "), join([@sprintf("%+6.2f", v) for v in hb], " "),
            mb(T, Tobs_g, okT .& dayt), mb(VP(T, RH), VPo, okV .& dayt), mb(VP(T, RH), VPo, okV .& nite),
            mb(T, Tobs_g, okT .& dayt .& EARLY), mb(T, Tobs_g, okT .& dayt .& LATE))
end
println("""
  Bias = model - measured (K; VP hPa). n54 posterior mean, full season from day 0, open loop except the CL row.
  split defaults = (lee + wind)/2 with vent_exp 1 (check: close to 'OL mean (current)').
  g = windward discharge multiplier, c = windward wind-pressure multiplier, e = opening curvature (both sides).""")
