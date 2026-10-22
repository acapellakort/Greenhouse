# =============================================================================
# crop.jl -- Aaron Velez's generic cucumber model (Marcelis 1994 in generic form) with his canopy/LAI extension
# =============================================================================
# Port of V2.2 generic_growth.jl (CUCUMBER_PARAMS_V2) and growth.jl `_grow_canopy!` with the n105 settings:
# node rate x M94 light factor, fruit sink held at stage 1 past maturity and picking by size, no temperature windows,
# no multipliers. States per m2 floor; one stem = one modelled plant (stem_density stems m-2).
# Daily order (Aaron, "Daily timing"): maintenance -> potentials of the start-of-day inventory -> allocation with
# reserve -> vegetative growth (leaf / non-leaf, area from the day's SLA) -> fruit development and node appearance
# (leaves, flowers, fruit set) -> harvest by size -> pruning to the retained-leaf target.

"M94 cucumber relations in generic form, from the four inferred physical parameters."
struct CucumberGeneric
    s::Float64; p::Float64; peak::Float64      # fruit sink shape
    A0::Float64; AT::Float64                   # sink amplitude A0 + AT T  (= AT (T - 10))
    kappa::Float64; base::Float64              # development clock kappa max(T - base, 0)
    beta0::Float64; beta1::Float64             # node appearance max(beta0 + beta1 T, 0)
    V20::Float64; VT::Float64                  # vegetative demand max(V20 + VT (T - 20), 0)
end

"""
    cucumber_generic(W_pot, E_fruit, beta1, V20, gp) -> CucumberGeneric

AT is chosen so that a fruit's total growth to maturity is W_pot for any E_fruit: AT = W_pot s (p - 1) / E_fruit.
"""
function cucumber_generic(W_pot, E, beta1, V20, gp)
    AT = W_pot * gp.sink_s * (gp.sink_p - 1) / E
    return CucumberGeneric(gp.sink_s, gp.sink_p, gp.sink_peak, -10.0 * AT, AT, 1.0 / E, gp.dev_base,
                           -beta1 * gp.app_base, beta1, V20, gp.VT)
end

@inline node_appearance(T, G::CucumberGeneric) = max(G.beta0 + G.beta1 * T, 0.0)
@inline function development_rate(T, G::CucumberGeneric)
    excess = T - G.base
    excess <= 0.0 && return 0.0
    return max(G.kappa * excess, 0.0)
end
@inline function fruit_sink(T, z, n, G::CucumberGeneric)
    (z < 0.0 || z > 1.0) && return 0.0
    amplitude = max(G.A0 + G.AT * T, 0.0)
    amplitude == 0.0 && return 0.0
    u = exp(-G.s * (z - G.peak)) / (G.p - 1.0)
    return max(Float64(n) * amplitude * u / (1.0 + u)^G.p, 0.0)
end
@inline vegetative_demand(T, G::CucumberGeneric) = max(G.V20 + G.VT * (T - 20.0), 0.0)

"M94 eq. 4 / Aaron's allocation with reserve: growth_i = P_i / P * min(supply + reserve, P)."
function allocate_daily(potentials, supply::Float64, reserve::Float64 = 0.0)
    demands   = collect(Float64, potentials)
    available = max(supply, 0.0) + max(reserve, 0.0)
    total     = sum(demands)
    total <= 0.0 && return zeros(length(demands)), available
    growth = demands .* (min(total, available) / total)
    remaining = available - sum(growth)
    remaining < 0.0 && (remaining = 0.0)
    return growth, remaining
end

"Seasonal specific leaf area [m2 g-1] (Heuvelink 1999 eq. 3); doy = day of year."
@inline seasonal_sla(doy, gp) = 1e-4 * (gp.sla_mean + gp.sla_amp * sin(2pi * mod((doy + gp.sla_phase) / 365.0, 1.0)))

mutable struct Fruit
    dev::Float64      # development stage (1 = maturity of M94)
    W::Float64        # cohort dry mass [g m-2]
    n::Float64        # fruits in the cohort [m-2]
    day_set::Int      # season day of fruit set
end

mutable struct CropState
    W_leaf::Float64; W_stem::Float64; W_root::Float64   # [g DM m-2]
    LAI::Float64
    leaves::Float64        # retained leaves per stem
    node::Float64          # node counter per stem (fractional progress kept)
    fruits::Vector{Fruit}
    reserve::Float64       # carbohydrate carried to the next day [g DM m-2]
    yield_DW::Float64      # cumulative harvest [g DM m-2]
    yield_FW::Float64      # cumulative harvest [kg FW m-2]
    rem_leaf_DW::Float64   # cumulative pruned leaf mass [g m-2]
    rem_LAI::Float64       # cumulative pruned leaf area [m2 m-2]
end

"Transplant state: LAI0, N0 nodes, leaf mass behind LAI0 at the transplant-day SLA, non-leaf mass in the vegetative split."
function init_crop_state(gp)
    Wl  = gp.LAI0 / seasonal_sla(gp.gen_doy0, gp)
    nl  = Wl * (1.0 - gp.gen_fL) / gp.gen_fL
    fst = gp.frac_stem / (gp.frac_stem + gp.frac_root)
    return CropState(Wl, fst * nl, (1.0 - fst) * nl, gp.LAI0, floor(gp.gen_N0), gp.gen_N0, Fruit[], 0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
    crop_step!(s, Pg, T, par_MJ, keep, day, gp, G) -> (n_harvested, sum of n x age of the harvested fruits)

One day: Pg [g CH2O m-2 d-1] gross assimilation, T [C] daily mean temperature, par_MJ [MJ m-2 d-1] daily canopy PAR
(NaN -> no light factor), keep = grower's kept-fruit share (NaN -> 1), day = season day (0 = transplant).
"""
function crop_step!(s::CropState, Pg::Float64, T::Float64, par_MJ::Float64, keep::Float64, day::Int, gp,
                    G::CucumberGeneric)
    fpar  = isnan(par_MJ) ? 1.0 : (1.0 - exp(-0.5 - 0.5 * par_MJ)) / (1.0 - exp(-0.5 - 0.5 * gp.par_ref))
    stems = gp.stem_density

    # 1. maintenance and supply (one construction cost asrq)
    W_fruit_tot = isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits)
    q10f = gp.Q10_resp ^ ((T - 25.0) / 10.0)
    Rm = q10f * (gp.k_m_leaf * s.W_leaf + gp.k_m_stem * s.W_stem + gp.k_m_root * s.W_root + gp.k_m_fruit * W_fruit_tot)
    dDM = max(Pg - Rm, 0.0) / gp.asrq

    # 2. potentials of the start-of-day inventory, allocation with reserve
    veg = vegetative_demand(T, G) * stems
    fsinks = [fruit_sink(T, min(f.dev, 1.0), f.n, G) for f in s.fruits]
    total_pot = veg + (isempty(fsinks) ? 0.0 : sum(fsinks))
    growths, s.reserve = allocate_daily(vcat([veg], fsinks), dDM, s.reserve)

    # 3. vegetative growth: leaf / non-leaf split; new area from the day's SLA on new leaf mass only
    vegDM = growths[1]
    leafG = gp.gen_fL * vegDM; nonleaf = vegDM - leafG
    fst = gp.frac_stem / (gp.frac_stem + gp.frac_root)
    s.W_leaf += leafG; s.W_stem += fst * nonleaf; s.W_root += (1.0 - fst) * nonleaf
    s.LAI    += seasonal_sla(gp.gen_doy0 + day, gp) * leafG
    for (i, f) in enumerate(s.fruits); f.W += growths[i + 1]; end

    # 4. development; node appearance (leaf at every node, flower from gen_flower_node); fruit set
    ddev = development_rate(T, G)
    for f in s.fruits; f.dev += ddev; end
    if day < gp.top_day
        old = s.node
        s.node += node_appearance(T, G) * fpar
        k0 = floor(Int, old) + 1; k1 = floor(Int, s.node)
        if k1 >= k0
            s.leaves += k1 - k0 + 1
            nfl = count(j -> j >= gp.gen_flower_node, k0:k1)
            if nfl > 0 && day >= gp.set_start_day
                set_frac = clamp((dDM / (total_pot + 1e-9) - gp.set_r_low) / (gp.set_r_high - gp.set_r_low), 0.0, 1.0)
                n_new = nfl * stems * set_frac * (isnan(keep) ? 1.0 : keep)
                n_new > 0.0 && push!(s.fruits, Fruit(0.0, 0.0, n_new, day))     # stage 0, grows from tomorrow
            end
        end
    end

    # 5. harvest by size (picking at harvest_FW grams fresh), safety cap at stage 3
    hw = gp.harvest_FW * gp.DMC_fruit
    keepf = Fruit[]; hn = 0.0; ha = 0.0
    for f in s.fruits
        if f.dev >= 3.0 || (f.n > 0.0 && f.W / f.n >= hw)
            s.yield_DW += f.W; s.yield_FW += (f.W / gp.DMC_fruit) / 1000.0
            hn += f.n; ha += f.n * (day - f.day_set)
        else
            push!(keepf, f)
        end
    end
    s.fruits = keepf

    # 6. pruning to the retained-leaf target (leaf mass and area removed in proportion)
    if s.leaves > gp.gen_leaf_target
        r = gp.gen_leaf_target / s.leaves
        s.rem_leaf_DW += (1.0 - r) * s.W_leaf; s.rem_LAI += (1.0 - r) * s.LAI
        s.W_leaf *= r; s.LAI *= r; s.leaves = gp.gen_leaf_target
    end
    return hn, ha
end
