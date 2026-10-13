# =============================================================================
# diagnose_sourcesink.jl  --  supply/demand regime diagnostic
# =============================================================================
# Runs the cucumber growth model at the POSTERIOR MEAN and records internal
# daily quantities that are not exposed by simulate_growth_measured:
#
#   ss          = dDM_pot / total_sink   (supply/demand ratio)
#                 < 1  → source-limited  (assimilate is the bottleneck)
#                 > 1  → sink-limited    (organs can't absorb all the assimilate)
#   dDM_pot     = net available dry matter per day [g DM/m2/d]
#   veg_sink    = vegetative demand [g DM/m2/d]
#   fruit_sink  = fruit demand [g DM/m2/d]
#   set_frac    = fraction of potential fruit set actually realised
#
# This tells us whether AiCU ever enters the source-limited branch where
# set_rate and Wf_max have *different* effects (not just their product).
#
# Run:  julia --project=invernaderosV2.2/GreenhouseSim \
#              invernaderosV2.2/GreenhouseSim/test/diagnose_sourcesink.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseSim.jl"))
using .GreenhouseSim
using CSV, DataFrames, JLD2, Plots, Statistics, Printf

# =============================================================================
# CONFIG — must match calibrate_growth.jl
# =============================================================================
INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
DATA         = joinpath(INV, "invernaderos2.0", "data")
TEAM         = "AiCU"
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(@__DIR__, "..", "configfiles", "cucumber_growth.json")
gh_csv       = joinpath(DATA, TEAM, "Greenhouse_climate.csv")
meteo_csv    = joinpath(DATA, "meteo.csv")
chain_path   = joinpath(@__DIR__, "posteriors", TEAM * "_growth_chain.jld2")
figdir       = joinpath(@__DIR__, "figures"); mkpath(figdir)

# =============================================================================
# load
# =============================================================================
params  = load_params(climate_json, crop_json)
gp_base = load_growth_params(growth_json)

println("Loading measured climate ($TEAM) ...")
season = load_measured_climate(gh_csv, meteo_csv, params)
println("  season length = ", season.ndays, " days")

println("Loading chain: ", chain_path)
JLD2.@load chain_path chain QoI burn_in chain_length

QoI_growth = ["node_rate", "veg_sink_max", "LAI_max", "set_rate", "Wf_max", "set_start_day"]
ng = length(QoI_growth)
@assert QoI[1:ng] == QoI_growth "QoI order mismatch"

post      = chain[burn_in:end, 1:ng]
post_mean = vec(mean(post, dims = 1))

println("\nPosterior mean parameters:")
for (n, v) in zip(QoI_growth, post_mean)
    @printf "  %-18s = %.4f\n" n v
end

# =============================================================================
# Diagnostic simulation: mirrors grow!() but records internals
# =============================================================================
# Replicate only what we need from growth.jl (engine is frozen; we don't touch it)
const _REF_DENSITY = 2.5
_sink_shape_d(dev) = 6.0 * dev * (1.0 - dev)

mutable struct _Fruit
    W   :: Float64
    dev :: Float64
    n   :: Float64    # number of fruits in cohort [per m2]
end

mutable struct _State
    W_leaf   :: Float64
    W_stem   :: Float64
    W_root   :: Float64
    LAI      :: Float64
    fruits   :: Vector{_Fruit}
    yield_DW :: Float64
    yield_FW :: Float64
    node     :: Float64
end

function _daily_assimilation(M::Matrix{Float64}, lai::Float64, p::NamedTuple, dt::Float64)
    size(M, 2) == 0 && return 0.0
    s = 0.0
    k_ext = 0.75   # fixed in JSON
    for j in 1:size(M, 2)
        s += GreenhouseSim.assimilation(M[1,j], M[2,j], M[3,j],
                                        (1.0 - exp(-k_ext * lai)) / k_ext, p)
    end
    return s * dt * 1e-3 * (30.0 / 44.0)
end

function _daily_mean_T(M::Matrix{Float64})
    size(M, 2) == 0 && return 20.0
    return sum(@view M[1,:]) / size(M, 2) - 273.15
end

function run_diagnostic(p, gp, season; LAI0 = 0.5)
    SLA   = gp.SLA
    W_leaf0 = LAI0 / SLA
    s = _State(W_leaf0, 8.0, 4.0, LAI0, _Fruit[], 0.0, 0.0, 0.0)

    days_out   = Int[]
    ss_out     = Float64[]
    dDM_pot_out = Float64[]
    veg_out    = Float64[]
    fruit_out  = Float64[]
    setfrac_out = Float64[]
    n_fruits_out = Float64[]

    dens    = gp.stem_density / _REF_DENSITY
    lai_cap = gp.LAI_max * dens

    for d in 0:(season.ndays - 1)
        M   = season.day_samples[d + 1]
        Pg  = _daily_assimilation(M, s.LAI, p, season.dt)
        Tc  = _daily_mean_T(M)
        dT  = max(Tc - gp.T_base, 0.0)

        # maintenance respiration
        W_fruit_tot = isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits)
        q10f = gp.Q10_resp ^ ((Tc - 25.0) / 10.0)
        Rm   = q10f * (gp.k_m_leaf * s.W_leaf + gp.k_m_stem * s.W_stem +
                       gp.k_m_root * s.W_root + gp.k_m_fruit * W_fruit_tot)

        dDM_pot = max(Pg - Rm, 0.0) / gp.asrq

        # development
        active = d < gp.top_day
        if active; s.node += gp.node_rate * dT; end
        ddev = dT / gp.DD_fruit
        for f in s.fruits; f.dev += ddev; end

        # sink strengths
        veg_sink = max(gp.veg_sink_max * dens * (1.0 - s.LAI / lai_cap), 0.0)
        fruit_sink_pre = isempty(s.fruits) ? 0.0 :
                         sum(f.n * gp.Wf_max * _sink_shape_d(f.dev) * ddev for f in s.fruits)
        total_sink_pre = veg_sink + fruit_sink_pre

        ss       = dDM_pot / (total_sink_pre + 1e-9)
        set_frac = active && d >= gp.set_start_day ?
                   clamp((ss - gp.set_r_low) / (gp.set_r_high - gp.set_r_low), 0.0, 1.0) : 0.0

        if active && d >= gp.set_start_day
            push!(s.fruits, _Fruit(0.0, 0.0, gp.set_rate * dens * set_frac))
        end

        # allocate
        fruit_sinks = [f.n * gp.Wf_max * _sink_shape_d(f.dev) * ddev for f in s.fruits]
        total_sink  = veg_sink + (isempty(fruit_sinks) ? 0.0 : sum(fruit_sinks))
        if total_sink > 0.0
            dDM = min(dDM_pot, total_sink)
            veg_DM = dDM * veg_sink / total_sink
            s.W_leaf += veg_DM * gp.frac_leaf
            s.W_stem += veg_DM * gp.frac_stem
            s.W_root += veg_DM * gp.frac_root
            for (i, f) in enumerate(s.fruits)
                f.W += dDM * fruit_sinks[i] / total_sink
            end
        end

        # harvest
        keep = _Fruit[]
        for f in s.fruits
            if f.dev >= 1.0
                s.yield_DW += f.W
                s.yield_FW += (f.W / gp.DMC_fruit) / 1000.0
            else
                push!(keep, f)
            end
        end
        s.fruits = keep

        # LAI update
        cap = min(lai_cap, gp.LAI_deleaf)
        if gp.SLA * s.W_leaf > cap; s.W_leaf = cap / gp.SLA; end
        s.LAI = gp.SLA * s.W_leaf

        push!(days_out, d)
        push!(ss_out, ss)
        push!(dDM_pot_out, dDM_pot)
        push!(veg_out, veg_sink)
        push!(fruit_out, fruit_sink_pre)
        push!(setfrac_out, set_frac)
        push!(n_fruits_out, Float64(length(s.fruits)))
    end

    return (; days = days_out, ss = ss_out, dDM_pot = dDM_pot_out,
              veg_sink = veg_out, fruit_sink = fruit_out,
              set_frac = setfrac_out, n_fruits = n_fruits_out,
              final_yield_FW = s.yield_FW)
end

# =============================================================================
# run at posterior mean
# =============================================================================
gp_fit = update_params(gp_base, QoI_growth, post_mean)
d = run_diagnostic(params, gp_fit, season; LAI0 = 0.5)

println("\nFinal yield (diagnostic sim): ", round(d.final_yield_FW, digits = 2), " kg/m2")

# fraction of season in each regime
src_lim = sum(d.ss .< 1.0) / season.ndays
snk_lim = sum(d.ss .>= 1.0) / season.ndays
println("\nRegime breakdown (by day):")
@printf "  source-limited (ss < 1) : %4.1f%%  → set_rate and Wf_max have DISTINCT effects here\n" 100*src_lim
@printf "  sink-limited   (ss ≥ 1) : %4.1f%%  → only set_rate × Wf_max matters here\n"            100*snk_lim
@printf "  median ss across season :  %.2f\n" median(d.ss)
@printf "  min / max ss            :  %.2f / %.2f\n" minimum(d.ss) maximum(d.ss)

# =============================================================================
# PLOTS
# =============================================================================

# --- Plot 1: supply/demand ratio over the season ---
p1 = plot(d.days, d.ss;
          lw = 2, lc = :steelblue, label = "ss = dDM_pot / total_sink",
          xlabel = "day of season", ylabel = "supply / demand ratio (ss)",
          title  = "Source–sink regime ($TEAM, posterior mean)",
          legend = :topright, ylims = (0, max(3.0, maximum(d.ss) * 1.1)))
hline!(p1, [1.0]; lw = 1.5, lc = :crimson, ls = :dash, label = "ss = 1 (regime boundary)")
hline!(p1, [gp_fit.set_r_low, gp_fit.set_r_high];
       lw = 1, lc = :darkorange, ls = :dot,
       label = "set_r_low / set_r_high (fruit-set thresholds)")
# shade source-limited region
src_days = filter(i -> d.ss[i] < 1.0, eachindex(d.days))
if !isempty(src_days)
    vspan!(p1, d.days[src_days[1]], d.days[src_days[end]];
           fillalpha = 0.08, lw = 0, fillcolor = :crimson,
           label = "source-limited zone")
end
savefig(p1, joinpath(figdir, "diag_sourcesink_ratio.png"))
println("\nSaved: diag_sourcesink_ratio.png")

# --- Plot 2: demand components stacked ---
p2 = plot(d.days, d.veg_sink;
          lw = 2, lc = :forestgreen, label = "veg sink",
          xlabel = "day of season", ylabel = "g DM m⁻² d⁻¹",
          title  = "Demand components vs supply ($TEAM, posterior mean)",
          legend = :topright)
plot!(p2, d.days, d.veg_sink .+ d.fruit_sink;
      lw = 2, lc = :darkorange, label = "veg + fruit sink (total)")
plot!(p2, d.days, d.dDM_pot;
      lw = 2, lc = :steelblue, ls = :dash, label = "dDM_pot (supply)")
savefig(p2, joinpath(figdir, "diag_sourcesink_components.png"))
println("Saved: diag_sourcesink_components.png")

# --- Plot 3: set_frac and fruit cohort count ---
p3 = plot(d.days, d.set_frac;
          lw = 2, lc = :purple, label = "set_frac",
          xlabel = "day of season", ylabel = "",
          title  = "Fruit-set fraction and cohort count ($TEAM)",
          legend = :topright)
plot!(twinx(p3), d.days, d.n_fruits;
      lw = 1.5, lc = :gray, ls = :dash, label = "cohorts on plant",
      ylabel = "# cohorts", legend = :topleft)
savefig(p3, joinpath(figdir, "diag_sourcesink_setfrac.png"))
println("Saved: diag_sourcesink_setfrac.png")

println("\nAll diagnostic figures in: ", figdir)
println("""
─────────────────────────────────────────────────────────
INTERPRETATION
  ss >> 1 throughout → crop always sink-limited
           → set_rate × Wf_max is the only identifiable combination
           → reparametrisation to (P=product, R=ratio) would show
             R_posterior ≈ R_prior (no information gained)

  ss crosses 1 repeatedly → both regimes present
           → R = set_rate / Wf_max carries information from
             source-limited periods; reparametrisation is worthwhile

  set_frac < 1 during source-limited periods
           → carbon stress visibly suppresses fruit set;
             Wf_max (demand per fruit) and set_rate have distinct
             effects on that suppression → parameters ARE separable
─────────────────────────────────────────────────────────
""")
