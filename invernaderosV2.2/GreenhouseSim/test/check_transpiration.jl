# =============================================================================
# check_transpiration.jl -- model transpiration vs irrigation - drain - crop water (Note 69; CHECK ONLY)
# =============================================================================
# Run from the GreenhouseSim folder (~1 min per model, forward only):
#     PRIORS=n54 julia --project=. test/check_transpiration.jl
#     PRIORS=n61 julia --project=. test/check_transpiration.jl
#
# Not part of any likelihood. The area the irrigation 'L per m2' refers to is not documented, so two bases are shown:
#   floor   (96 m2):  water per floor m2 = logged value
#   growing (76.8 m2): water per floor m2 = logged value x 76.8/96
# Observed transpiration per floor m2 = net water - harvested fruit water - vegetative growth water:
#   fruit water = daily harvest (kg FW/m2 of production area) x 0.96 x 76.8/96 (production assumed uniform over gutters)
#   vegetative water = 0.06 L/m2 growing area/day (estimate, Note 68)
# Model transpiration = sum of p1 = q1 (Psat(T1) - V1) over the day (climate_rhs, stomata_fix form), per floor m2.
# Weekly sums remove most of the slab-storage swing (~ +-0.5 L/m2).
# =============================================================================

ENV["CL_MODE"] = "setup"; haskey(ENV, "PRIORS") || (ENV["PRIORS"] = "n61")
haskey(ENV, "RELTOL") || (ENV["RELTOL"] = "1e-4")
include(joinpath(@__DIR__, "run_climate_calibration_cl.jl"))

parts = Vector{Matrix{Float64}}()
for s in (1, 2)
    f = joinpath(POSTDIR, "AiCU_climcal_cl_blocks_$(PRIORS)_seed$(s).jld2"); fc = replace(f, ".jld2" => "_ckpt.jld2")
    src = isfile(f) ? f : isfile(fc) ? fc : ""; isempty(src) && continue
    ch = JLD2.load(src)["chain"]; size(ch, 1) > 31000 || continue
    push!(parts, ch[30001:end, 1:NM]); println("seed $s: ", basename(src))
end
isempty(parts) && error("no $(PRIORS) chains in $POSTDIR")
pmean = vec(mean(reduce(vcat, parts), dims = 1)); pT = update_params(params_cl, CQ_MODEL, pmean)
# Note 71: second run with the crop physiology seeing sun-through-screens + lamp PAR (phys_light) and the
# leaf-area-scaled photosynthesis (photo_lai); same posterior-mean parameters (not re-fitted)
pL = update_params(pT, ["phys_light", "photo_lai", "lamp_par_W"], [1.0, 1.0, 187.0 / 4.6])

function p1_flux(T1, T2, C1, V1, I9, LAI, p)        # kg m-2 s-1, as climate_rhs with stomata_fix = 1
    (; gamma, gamma2, gamma3, gamma4, gamma5, delta1, delta2, delta3, delta4, delta5, delta6, delta7, rho3, alpha5) = p
    q2  = Pws(T1)
    q7  = 1 / (1 + exp((I9 - delta1) / gamma5))
    q10 = (I9 + delta2) / (I9 + delta3)
    q4  = 1 + (delta5 + (delta4 - delta5) * q7) * (C1 * T2 / 536.4 - 200)^2
    q5  = 1 + (delta7 + (delta6 - delta7) * q7) * (q2 - V1)^2
    rs  = gamma4 * q10 * q4 * q5
    q1  = 2 * rho3 * alpha5 * LAI / (gamma * gamma2 * (gamma3 + rs))
    return q1 * (q2 - V1)
end
trans9 = (1 - pT.eta1) * pT.tau1 * pT.eta2
Ig = A(Ig_v).(predtimes); LAIt = A(LAI_row).(predtimes)
BSt = A(Blackv).(predtimes); ESt = A(U1v).(predtimes); LMt = A(U12v).(predtimes)
function emod(pX, phys)
    us_ = simdays(pX, days, u0_cl, RELTOL); us_ === nothing && error("solve failed")
    E = zeros(N_DAYS)
    for d in 0:(N_DAYS - 1), j in 0:(NPD - 2)
        c = col_of(d, j)
        Ipar = phys ? trans9 * Ig[c] * (1 - 0.98 * clamp(BSt[c], 0, 1)) * (1 - 0.25 * clamp(ESt[c], 0, 1)) + pX.lamp_par_W * clamp(LMt[c], 0, 1) :
                      trans9 * Ig[c]
        E[d + 1] += p1_flux(us_[1, c], us_[2, c], us_[3, c], us_[4, c], Ipar, LAIt[c], pX) * 300.0
    end
    E
end
Emod  = emod(update_params(pT, ["phys_light", "photo_lai"], [0.0, 0.0]), false)   # sun only (pre-Note 72 behaviour)
EmodL = emod(pL, true)

irr = DataFrame(CSV.File(joinpath(DATA, TEAM, "Irrigation.csv")))
num(x) = x isa Number ? Float64(x) : NaN
net = fill(NaN, N_DAYS)
for r in eachrow(irr)
    d = Int(floor(num(r.time) - DAY0_SERIAL)); 0 <= d < N_DAYS || continue
    net[d + 1] = num(r.water) - (isnan(num(r.drain)) ? 0.0 : num(r.drain))
end
harv = [cum_at(days_obs, y_obs, d) - cum_at(days_obs, y_obs, d - 1) for d in 0:(N_DAYS - 1)]; harv[1] = 0.0
fruitw = harv .* 0.96 .* (76.8 / 96)
vegw   = fill(0.06 * 76.8 / 96, N_DAYS)
Eobs(basis) = net .* (basis / 96) .- fruitw .- vegw

println("\n", "="^104)
@printf("  %s transpiration check (posterior mean, %s); L per floor m2 per week\n", PRIORS, OLP ? "open loop" : "closed loop")
@printf("  %-10s | %8s %8s | %8s %8s | %8s %8s | %6s %6s | %6s %6s\n", "days", "net(96)", "fruit", "Eobs 96", "Eobs 76.8",
        "Emod", "Emod+L", "r96", "r77", "r96+L", "r77+L")
W = Tuple{Int,Int}[(a, min(a + 6, N_DAYS - 1)) for a in 0:7:(N_DAYS - 1)]
wkd = Dict{Int,Vector{Float64}}()
for (a, b) in W
    ix = (a:b) .+ 1; ok = .!isnan.(net[ix]); count(ok) >= 6 || continue
    e96 = sum(Eobs(96.0)[ix][ok]); e77 = sum(Eobs(76.8)[ix][ok]); em = sum(Emod[ix][ok]); eL = sum(EmodL[ix][ok])
    wkd[a] = [e96, e77, em, eL]
    @printf("  d%-3d-%-4d | %8.2f %8.2f | %8.2f %8.2f | %8.2f %8.2f | %6.2f %6.2f | %6.2f %6.2f\n", a, b, sum(net[ix][ok]),
            sum(fruitw[ix][ok]), e96, e77, em, eL, em / e96, em / e77, eL / e96, eL / e77)
end
println("  ", "-"^102)
for (lo, hi) in [(0, 8), (9, 48), (49, 100), (101, 114), (0, 114)]
    ks = [k for k in keys(wkd) if lo <= k <= hi]; isempty(ks) && continue
    s = reduce(+, [wkd[k] for k in ks])
    @printf("  d%-3d-%-4d | model/obs  sun only: %.2f floor, %.2f growing | sun+lamps: %.2f floor, %.2f growing  [%d wk]\n",
            lo, hi, s[3] / s[1], s[3] / s[2], s[4] / s[1], s[4] / s[2], length(ks))
end
ks = sort(collect(keys(wkd)))
@printf("  weekly correlation (model, observed): sun only %.2f, sun+lamps %.2f\n", cor([wkd[k][3] for k in ks], [wkd[k][1] for k in ks]),
        cor([wkd[k][4] for k in ks], [wkd[k][1] for k in ks]))
println("""
  Read: ratio ~1 on one basis = the model transpiration is consistent with the water use on that area basis.
  Model LAI = interpolated derived LAI (2.5 plateau, per growing area) with no 0.8 floor factor. A ratio that
  drifts with the season points at LAI or stomata; a constant ratio at the area basis or the fruit/veg terms.""")
