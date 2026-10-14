# =============================================================================
# diagnose_lamp_screen.jl -- does the missing lamp heat explain the day offset (AiCU) and the late-season cold
#                           nights (other teams)? (Note 54; forward only)
# =============================================================================
# Run from the GreenhouseSim folder (forward runs only, ~2 min per team, no MCMC):
#     for T in AiCU Croperators DeepGreens iGrow Sonoma "Reference(Growers)"; do
#         TEAM="$T" julia --project=. test/diagnose_lamp_screen.jl; done
#     for T in Croperators DeepGreens Sonoma; do                       # second heating circuit as well
#         PIPE2=1 TEAM="$T" julia --project=. test/diagnose_lamp_screen.jl; done
#
# Physics = AiCU n49b pooled posterior. Four runs per team: {CL, OL} x {V2.0 lamp heat (29+4 W m-2), lamp_full
# (117 W m-2 electrical; 45 % x (1-exp(-0.7 LAI)) to the canopy, rest to the air)}. Extra columns split the
# T bias by lamp state (AssimLight > 50 %): if the missing lamp heat is the cause, the cold bias sits in the
# "lamp on" columns under V2.0 and shrinks there under lamp_full, while the "lamp off" columns barely move.
# Caveat: n49b was calibrated with the V2.0 lamp heat, so its other parameters may already absorb part of it.
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n49b")   # Note 55: PRIORS=n54 uses the n54 posterior
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2")
    isfile(f) || error("missing $f")
    ch = JLD2.load(f)["chain"]; push!(parts, ch[30001:end, 1:NM])
end
pm = vec(mean(reduce(vcat, parts), dims = 1))
pP = update_params(params_cl, CQ_MODEL, pm)

okT = .!isnan.(Tobs_g); okT[1] = false; okV = okT .& .!isnan.(RHobs_g); okU = .!isnan.(U8obs); okU[1] = false
inset(ds) = in.(cd_, Ref(Set(ds)))
ALLD = inset(collect(0:114))
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
Igc = A(Ig_v).(predtimes)
Toc = A(Tout_v).(predtimes)

m = okT .& dayt
@printf("\nTEAM %s: measured day T - T_out = %+.2f K, measured night T - T_out = %+.2f K, mean logged U8 day %.3f / night %.3f\n",
        TEAM, mean(Tobs_g[m] .- Toc[m]), mean(Tobs_g[okT .& nite] .- Toc[okT .& nite]),
        mean(U8obs[okU .& dayt]), mean(U8obs[okU .& nite]))

U12c = A(U12v).(predtimes); lon = U12c .> 0.5
@printf("PIPE2 = %s (PIPE2_R %.2f); lamp on: %.0f %% of day, %.0f %% of night samples\n", PIPE2, PIPE2_R, 100mean(lon[dayt]), 100mean(lon[nite]))
println("\n", "="^200)
@printf("  %-24s | %6s | %s | %s | %6s %6s %6s %7s | %s\n", "TEAM / mode", "hRMS",
        join([@sprintf("%6s", "n$(a)") for (a, b) in PER], " "), join([@sprintf("%6s", "d$(a)-$(b)") for (a, b) in HB_], " "),
        "T day", "VPday", "VPnit", "U8 day", "sensor test a + b*I/100   | day L-on L-off | nit L-on L-off")
println("  ", "-"^198)
# n49b posterior: base = V2.0 lamp heat, variant = lamp_full. n54 posterior (lamp_full fitted on): base = V2.0, variant = n54 as fitted.
p0 = update_params(pP, ["lamp_full"], [0.0]); pL = update_params(pP, ["lamp_full"], [1.0])
@printf("posterior: %s (pooled seeds 1-2, steps 30001+)\n", PRIORS)
for (mode, p) in (("CL", p0), ("CL lamp_full", pL), ("OL", update_params(p0, ["vent_ctrl"], [0.0])),
                  ("OL lamp_full", update_params(pL, ["vent_ctrl"], [0.0])))
    us = simdays(p, days, u0_cl, RELTOL)
    us === nothing && (println("  $TEAM $mode: solver failed"); continue)
    T = us[2, :]; RH = rhf.(us[2, :], us[4, :]); U = startswith(mode, "CL") ? clamp.(us[7, :], 0.0, 1.0) : U8obs
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    r = Tobs_g[m] .- T[m]; ab = hcat(ones(count(m)), Igc[m] ./ 100) \ r
    lb = [mb(T, Tobs_g, okT .& w .& l) for (w, l) in ((dayt, lon), (dayt, .!lon), (nite, lon), (nite, .!lon))]
    @printf("  %-24s | %6.3f | %s | %s | %+6.2f %+6.2f %+6.2f %+7.3f | a %+.2f  b %+.3f | %+6.2f %+6.2f | %+6.2f %+6.2f\n", "$(TEAM) $(mode)",
            hourly_rmse(T, Tobs_g, okT .& ALLD), join([@sprintf("%+6.2f", v) for v in nb], " "),
            join([@sprintf("%+6.2f", v) for v in hb], " "), mb(T, Tobs_g, okT .& ALLD .& dayt),
            mb(VP(T, RH), VPo, okV .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& nite),
            startswith(mode, "CL") ? mb(U, U8obs, okU .& ALLD .& dayt) : 0.0, ab[1], ab[2], lb...)
end
println("""
  Bias = model - measured (K). L-on / L-off: samples with lamps on / off. V2.0 lamp heat = 33 W m-2 vs 117 W m-2
  electrical. Read: does lamp_full remove the late-season night cold (n49..n110) of the lamp-heavy teams and the
  AiCU day deficit in the L-on column, without making the L-off columns move?""")
