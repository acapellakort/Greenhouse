# =============================================================================
# diagnose_teams_day.jl -- is the ~1 K daytime offset specific to AiCU, or common to all AGC-2018 compartments? (Note 54)
# =============================================================================
# Run from the GreenhouseSim folder, once per team (forward runs only, ~1 min each, no MCMC):
#     for T in AiCU Croperators DeepGreens iGrow Sonoma "Reference(Growers)"; do
#         TEAM="$T" julia --project=. test/diagnose_teams_day.jl; done
#
# Physics = the AiCU n49b pooled posterior (Note 53), applied unchanged to each team's own inputs (weather, screens,
# pipes, setpoints, LAI from its CropManagement.csv). Two modes:
#   CL  closed-loop vents with the AiCU controller parameters (grower-specific: a mismatch is expected)
#   OL  open loop: the team's LOGGED vent opening is replayed (vent_ctrl = 0) -> pure energy balance, no controller
# The OL day bias by hour of day is the key number: if every compartment shows the same ~1 K afternoon cold bias
# with its own logged vents, the offset is common to the measurement set-up or to missing physics shared by all
# compartments; if only AiCU shows it, it is specific to that compartment (its sensor or its management).
# =============================================================================

ENV["CL_MODE"] = "setup"; ENV["PRIORS"] = "n49b"
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_n49b_seed$(s).jld2")
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

println("\n", "="^170)
@printf("  %-14s | %6s | %s | %s | %6s %6s %6s %7s | %s\n", "TEAM / mode", "hRMS",
        join([@sprintf("%6s", "n$(a)") for (a, b) in PER], " "), join([@sprintf("%6s", "d$(a)-$(b)") for (a, b) in HB_], " "),
        "T day", "VPday", "VPnit", "U8 day", "sensor test a + b*I/100")
println("  ", "-"^168)
for (mode, p) in (("CL", pP), ("OL", update_params(pP, ["vent_ctrl"], [0.0])))
    us = simdays(p, days, u0_cl, RELTOL)
    us === nothing && (println("  $TEAM $mode: solver failed"); continue)
    T = us[2, :]; RH = rhf.(us[2, :], us[4, :]); U = mode == "CL" ? clamp.(us[7, :], 0.0, 1.0) : U8obs
    nb = [mb(T, Tobs_g, okT .& nite .& inset(collect(a:b))) for (a, b) in PER]
    hb = [mb(T, Tobs_g, okT .& dayt .& (hrc .>= a) .& (hrc .< b)) for (a, b) in HB_]
    r = Tobs_g[m] .- T[m]; ab = hcat(ones(count(m)), Igc[m] ./ 100) \ r
    @printf("  %-14s | %6.3f | %s | %s | %+6.2f %+6.2f %+6.2f %+7.3f | a %+.2f  b %+.3f\n", "$(TEAM) $(mode)",
            hourly_rmse(T, Tobs_g, okT .& ALLD), join([@sprintf("%+6.2f", v) for v in nb], " "),
            join([@sprintf("%+6.2f", v) for v in hb], " "), mb(T, Tobs_g, okT .& ALLD .& dayt),
            mb(VP(T, RH), VPo, okV .& ALLD .& dayt), mb(VP(T, RH), VPo, okV .& ALLD .& nite),
            mode == "CL" ? mb(U, U8obs, okU .& ALLD .& dayt) : 0.0, ab[1], ab[2])
end
println("""
  OL = logged vents replayed (no controller): the cleanest test of the energy balance.
  Compare the OL day-by-hour row across teams. Same physics everywhere; LAI, screens, pipes and vents are each
  team's own logs. CL rows for other teams use AiCU's controller parameters and are expected to be worse.""")
