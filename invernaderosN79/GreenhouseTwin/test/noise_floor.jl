# =============================================================================
# test/noise_floor.jl  --  numerical noise diagnostic for n79 GreenhouseTwin
# =============================================================================
# Reproduces the Note 109 noise table and isolates the contribution of the
# 5-min input discontinuities by adding tstops to the ODE solver.
# No production code is modified. A local simdays_tstops is defined here.
#
# Run from  invernaderosN79/GreenhouseTwin:
#   CROPTAG=n105 julia --project=. test/noise_floor.jl 2>&1 | tee noise_floor.log
#
# Metric: perturb PERT_PARAM by k * 1e-14 (k = 1..N_PERT), then measure:
#   (a) max|ΔT2| -- maximum absolute T2 difference across all TRAIN block days [K]
#   (b) energy sd and range over the N_PERT + 1 evaluations (baseline included)
# =============================================================================

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "GreenhouseTwin.jl"))
using .GreenhouseTwin
using OrdinaryDiffEq          # needed for simdays_tstops (ODEProblem, Tsit5, solve)
using Printf, Statistics

# ── identical setup to calibrate_climate.jl ──────────────────────────────────
const GT      = GreenhouseTwin
const INV     = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CJ      = joinpath(@__DIR__, "..", "configfiles", "climate_n79.json")
const CROPTAG = get(ENV, "CROPTAG", "n105")
const CROPCSV = let
    n = joinpath(@__DIR__, "posteriors", "AiCU_crop_series_$CROPTAG.csv")
    o = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "test", "posteriors",
                 "AiCU_crop_series_$CROPTAG.csv")
    isfile(n) ? n : o
end

const p0      = load_params(CJ)
const PR      = priors_of(CJ)
const Q       = vcat(GT.CLIMATE_QOI, GT.CLIMATE_NOISE)
const S       = load_agc_season(INV; crop_series = CROPCSV)
const BL      = climate_blocks(S)
const centre  = [PR[q][2] for q in Q]
const p_ctr   = update_params(p0, GT.CLIMATE_QOI, centre[1:length(GT.CLIMATE_QOI)])

const PERT_PARAM = "nu1"   # ventilation discharge: active day and night, not inside a switch
const EPS        = 1e-14   # relative size of each perturbation
const N_PERT     = 10      # number of perturbations

# ── Variant 1: simdays with tstops at every 5-min input-step boundary ────────
# Identical to GT.simdays; the only change is `tstops` passed to `solve`.
# Each input function is a piecewise constant step function on the 300-s grid;
# stopping the solver at those boundaries prevents steps that straddle a jump.
function simdays_tstops(p, season, days, u0, reltol;
                        tspan = (0.0, 86400.0),
                        teval = 0.0:300.0:86400.0)
    ts  = collect(0.0:300.0:86400.0)   # 289 stop points = every 5-min input boundary
    npd = length(teval)
    us  = Matrix{Float64}(undef, length(u0), 1 + length(days) * npd)
    us[:, 1] = u0;  u = copy(u0);  col = 2
    for d in days
        ctx = SimContext(p, season.weather, season.controls,
                         season.t0 + 86400.0 * d, 86400.0 * d)
        sol = try
            solve(ODEProblem(rhs!, u, tspan, ctx), Tsit5();
                  saveat = collect(teval), tstops = ts, reltol = reltol)
        catch
            return nothing
        end
        length(sol.t) == npd || return nothing
        A = Array(sol);  us[:, col:col+npd-1] = A;  col += npd;  u = A[:, end]
    end
    return us
end

# ── Pluggable energy and T2-noise functions ───────────────────────────────────
# These mirror block_sse / climate_energy from climate_calibration.jl but
# accept any simdays-like function as an argument.

function _bsse(v, bins, obs)     # bin-mean squared error
    s = 0.0
    for (ix, o) in zip(bins, obs)
        m = 0.0; @inbounds for i in ix; m += v[i]; end
        d = m / length(ix) - o;  s += d * d
    end
    return s
end

function block_sse_v(p, b, rtol, sfn)
    us = sfn(p, S, b.days, b.u0, rtol)
    us === nothing && return (Inf, Inf, Inf)
    T2 = us[2, :];  C1 = us[3, :];  RH = rhf.(T2, us[4, :])
    sT = _bsse(T2, b.binsT, b.obsT)
    sR = _bsse(RH, b.binsR, b.obsR)
    sC = _bsse(C1, b.binsC, b.obsC)
    return (isfinite(sT) ? sT : Inf, isfinite(sR) ? sR : Inf, isfinite(sC) ? sC : Inf)
end

function energy_v(z, sfn; rtol)
    nm = length(GT.CLIMATE_QOI)
    p  = update_params(p0, GT.CLIMATE_QOI, z[1:nm])
    sT = sR = sC = 0.0
    for b in BL
        a, r, c = block_sse_v(p, b, rtol, sfn)
        (isfinite(a) && isfinite(r) && isfinite(c)) || return Inf
        sT += a;  sR += r;  sC += c
    end
    nT = sum(length(b.obsT) for b in BL)
    nR = sum(length(b.obsR) for b in BL)
    nC = sum(length(b.obsC) for b in BL)
    σT, σR, σC = z[nm+1], z[nm+2], z[nm+3]
    return nT*log(σT) + sT/(2σT^2) + nR*log(σR) + sR/(2σR^2) + nC*log(σC) + sC/(2σC^2)
end

# ── Main measurement ──────────────────────────────────────────────────────────
pidx = findfirst(==(PERT_PARAM), Q)
pidx === nothing   && error("$PERT_PARAM not found in Q")
pidx > length(GT.CLIMATE_QOI) && error("$PERT_PARAM is a noise term; choose a QOI parameter")
const nm = length(GT.CLIMATE_QOI)

println("="^72)
println("noise_floor.jl  —  n79 GreenhouseTwin numerical noise diagnostic")
@printf("  perturbing  %-14s  ε = %.0e × k  (k = 0..%d)\n", PERT_PARAM, EPS, N_PERT)
@printf("  crop series: %s\n", basename(CROPCSV))
println("="^72)

println("\nWarm-up (1 energy eval at reltol 1e-3 for each variant) …")
energy_v(centre, simdays;        rtol = 1e-3)
energy_v(centre, simdays_tstops; rtol = 1e-3)
println("done.")

struct NStat; t2_max::Float64; e_sd::Float64; e_rng::Float64; end

function run_variant(label, rtol, sfn)
    @printf("\n[%s  reltol %.0e]\n", label, rtol)
    # pre-compute baseline trajectories once (avoid re-running for each perturbation)
    ub_all = [sfn(p_ctr, S, b.days, b.u0, rtol) for b in BL]
    any(isnothing, ub_all) && error("baseline solve failed for variant '$label'")
    e0 = energy_v(centre, sfn; rtol)
    @printf("  k= 0  energy = %11.3f  (baseline)\n", e0)
    t2s = Float64[];  ens = Float64[e0]
    for k in 1:N_PERT
        z      = copy(centre);  z[pidx] *= (1.0 + k * EPS)
        p_pert = update_params(p0, GT.CLIMATE_QOI, z[1:nm])
        # T2 noise: max |ΔT2| across all TRAIN blocks
        maxd = 0.0
        for (ub, b) in zip(ub_all, BL)
            up = sfn(p_pert, S, b.days, b.u0, rtol)
            if up === nothing; maxd = Inf; break; end
            maxd = max(maxd, maximum(abs.(up[2, :] .- ub[2, :])))
        end
        en = energy_v(z, sfn; rtol)
        push!(t2s, maxd);  push!(ens, en)
        @printf("  k=%2d  |ΔT2|_max = %7.4f K   energy = %11.3f   Δ = %+.3f\n",
                k, maxd, en, en - e0)
    end
    NStat(maximum(t2s), std(ens), maximum(ens) - minimum(ens))
end

variants = [
    ("0:  standard (no tstops)", 1e-4, simdays),
    ("1:  + tstops at 5-min   ", 1e-4, simdays_tstops),
    ("0b: standard (no tstops)", 1e-6, simdays),
    ("1b: + tstops at 5-min   ", 1e-6, simdays_tstops),
]

stats = NStat[]
for (lbl, rtol, sfn) in variants
    t0 = time()
    push!(stats, run_variant(lbl, rtol, sfn))
    @printf("  → %.1f s\n", time() - t0)
end

println("\n" * "="^72)
println("SUMMARY TABLE")
println("="^72)
@printf("%-30s  %8s  %13s  %10s  %13s\n",
        "variant", "reltol", "max|ΔT2| (K)", "energy sd", "energy range")
for ((lbl, rtol, _), st) in zip(variants, stats)
    @printf("%-30s  %8.0e  %13.4f  %10.2f  %13.2f\n",
            lbl, rtol, st.t2_max, st.e_sd, st.e_rng)
end

println()
println("Reference (Note 109, Mac, reltol 1e-4, 1e-14 perturbation of one parameter):")
println("  max|ΔT2| up to 1.7 K,  energy sd 13-27,  energy range 37-62")
println("Reference (reltol 1e-6):")
println("  energy sd 0.36-0.64,  energy range 1.0-2.1")
