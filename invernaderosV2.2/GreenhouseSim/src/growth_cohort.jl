# =============================================================================
# growth_cohort.jl -- leaf-age cohort canopy + node-driven flowering (crop v2)
# =============================================================================
# Fixes the big-leaf model's inability to represent de-leafing:
#   * leaves are AGE COHORTS produced at the head (node rate) -> no regrowth churn;
#   * de-leafing removes the OLDEST (bottom, shaded) cohorts, shedding their
#     maintenance respiration while total interception (saturated) barely drops;
#   * flowering is NODE-DRIVEN: each node's flowers set (carbon-regulated) into
#     fruit cohorts carrying their own degree-day age -> growth and harvest.
# Photosynthesis stays big-leaf on TOTAL LAI (fast, keeps the climate CO2 coupling);
# the cohorts carry age, respiration, and de-leafing. Reuses Fruit, _sink_shape,
# REF_DENSITY, lai_eff, daily_assimilation_measured from growth.jl / measured.jl.

mutable struct LeafCohort
    area::Float64   # LAI contribution [m2 leaf / m2 ground]
    age::Float64    # degree-days since appearance
end

mutable struct CohortState
    leaves::Vector{LeafCohort}
    W_stem::Float64
    W_root::Float64
    fruits::Vector{Fruit}
    yield_DW::Float64
    yield_FW::Float64
    node::Float64
    LAI::Float64
end

leaf_mass(s::CohortState, SLA) = s.LAI / SLA         # total leaf DM [g/m2]

"""
    canopy_light_profile(leaves, k_ext) -> (order, L_mid)

Sort cohorts young(top)->old(bottom) and return, for each in that order, the
cumulative-LAI depth to the cohort mid-point. Light reaching a cohort is
I0*exp(-k_ext*L_mid): the vertical light gradient the big-leaf model averages away.
"""
function canopy_light_profile(leaves::Vector{LeafCohort}, k_ext::Float64)
    order = sortperm(leaves, by = c -> c.age)      # youngest (top) first
    Lmid  = Vector{Float64}(undef, length(order))
    Lup   = 0.0
    @inbounds for (r, i) in enumerate(order)
        a = leaves[i].area
        Lmid[r] = Lup + 0.5 * a
        Lup    += a
    end
    return order, Lmid
end

"""
    daily_assimilation_layered(M, leaves, p, dt, k_ext) -> g CH2O m^-2 d^-1

Layered canopy photosynthesis: integrate the FvCB LEAF rate down the light
gradient (midpoint rule over the age-sorted cohorts) instead of running it at
top light and scaling by lai_eff. Reduces to the big-leaf value when the leaf
rate is linear in light (light-limited crop); diverges (correctly lower) once
top leaves light-saturate, and it charges deep, shaded leaves ~zero fixation.
"""
function daily_assimilation_layered(M::Matrix{Float64}, leaves::Vector{LeafCohort},
                                    p::NamedTuple, dt::Float64, k_ext::Float64 = K_EXT)
    (isempty(leaves) || size(M, 2) == 0) && return 0.0
    order, Lmid = canopy_light_profile(leaves, k_ext)
    tot = 0.0
    @inbounds for j in 1:size(M, 2)
        T = M[1, j]; C = M[2, j]; Iw = M[3, j]
        VPD = size(M, 1) >= 4 ? M[4, j] : 0.0     # VPD row optional (0 -> constant-stomata path)
        for (r, i) in enumerate(order)
            a = leaves[i].area
            a <= 0.0 && continue
            Iloc = Iw * exp(-k_ext * Lmid[r])          # light reaching this cohort
            tot += assimilation(T, C, Iloc, a, p; VPD = VPD)   # per-cohort leaf area, local light
        end
    end
    return tot * dt * 1e-3 * (30.0 / 44.0)             # -> g CH2O m^-2 d^-1
end

function init_cohort_state(gp; LAI0 = 0.5, stem0 = 8.0, root0 = 4.0)
    return CohortState([LeafCohort(LAI0, 0.0)], stem0, root0, Fruit[], 0.0, 0.0, 0.0, LAI0)
end

"""
    grow_cohort!(s, Pg, Tmean_C, gp, day; deleaf_target=Inf, generic_params=nothing)

Advance one day. `Pg` = daily gross assimilate [g CH2O/m2] (big-leaf, total LAI).
`deleaf_target` = working LAI the grower holds by removing the oldest leaves.

`generic_params`: pass a `GenericParameters` instance (e.g. `CUCUMBER_PARAMS`)
to use the Marcelis/Heuvelink Bell-curve sink and affine appearance rate in place
of the legacy smoothstep sink and capacity-fill vegetative demand. Default
`nothing` preserves the original V2.2 behaviour exactly.
"""
function grow_cohort!(s::CohortState, Pg::Float64, Tmean_C::Float64, gp, day::Int;
                      deleaf_target = Inf, dens = gp.stem_density / REF_DENSITY,
                      q10_resp = NaN,
                      generic_params::Union{Nothing, GenericParameters} = nothing)
    dT   = max(Tmean_C - gp.T_base, 0.0)
    q10f = isnan(q10_resp) ? gp.Q10_resp ^ ((Tmean_C - 25.0) / 10.0) : q10_resp  # diurnal-corrected if provided

    # --- maintenance respiration on all standing tissue -----------------------
    Wl = leaf_mass(s, gp.SLA)
    Wf = isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits)
    Rm = q10f * (gp.k_m_leaf * Wl + gp.k_m_stem * s.W_stem + gp.k_m_root * s.W_root + gp.k_m_fruit * Wf)
    dDM_pot = max(Pg - Rm, 0.0) / gp.asrq

    # --- development: nodes appear (head), cohorts & fruits age ----------------
    active = day < gp.top_day
    dnode  = active ? gp.node_rate * dT : 0.0
    s.node += dnode

    # Fruit development increment: generic clock or legacy linear thermal time
    if isnothing(generic_params)
        ddev = dT / gp.DD_fruit
    else
        ddev = generic_development_rate(Tmean_C, generic_params)
    end
    for f in s.fruits; f.dev += ddev; end
    for c in s.leaves; c.age += dT;  end

    if isnothing(generic_params)
        # --- Legacy path (V2.2 smoothstep) ---

        # sinks: node-driven LEAF demand (no capacity-fill churn) + fruit
        leaf_area_demand = gp.leaf_per_node * dens * dnode          # new LAI wanted at the head
        veg_sink = (leaf_area_demand / gp.SLA) / max(gp.frac_leaf, 1e-6)   # total veg DM (leaf is frac_leaf share)
        fruit_sink_pre = isempty(s.fruits) ? 0.0 :
                         sum(f.n * gp.Wf_max * _sink_shape(f.dev) * ddev for f in s.fruits)

        # node-driven flowering: flowers per node set (carbon-regulated)
        if active && day >= gp.set_start_day && dnode > 0.0
            ss = dDM_pot / (veg_sink + fruit_sink_pre + 1e-9)
            set_frac = clamp((ss - gp.set_r_low) / (gp.set_r_high - gp.set_r_low), 0.0, 1.0)
            n_set = gp.set_rate * dens * dnode * set_frac            # fruits set at this node flush
            n_set > 0.0 && push!(s.fruits, Fruit(0.0, 0.0, n_set))
        end

        # allocate assimilate: source- or sink-limited
        fruit_sinks = [f.n * gp.Wf_max * _sink_shape(f.dev) * ddev for f in s.fruits]
        total_sink  = veg_sink + (isempty(fruit_sinks) ? 0.0 : sum(fruit_sinks))
        if total_sink > 0.0
            dDM    = min(dDM_pot, total_sink)
            veg_DM = dDM * veg_sink / total_sink
            new_leaf_area = veg_DM * gp.frac_leaf * gp.SLA           # new top cohort
            new_leaf_area > 0.0 && push!(s.leaves, LeafCohort(new_leaf_area, 0.0))
            s.W_stem += veg_DM * gp.frac_stem
            s.W_root += veg_DM * gp.frac_root
            for (i, f) in enumerate(s.fruits); f.W += dDM * fruit_sinks[i] / total_sink; end
        end

    else
        # --- Generic path (Marcelis 1994 / Heuvelink 1996 Bell sigmoid) ---

        # Node-driven leaf demand (same node rate as legacy, generic sink replaces fruit shape)
        leaf_area_demand = gp.leaf_per_node * dens * dnode
        veg_sink = (leaf_area_demand / gp.SLA) / max(gp.frac_leaf, 1e-6)

        fruit_sink_pre = isempty(s.fruits) ? 0.0 :
                         sum(generic_sink(Tmean_C, f.dev, f.n, generic_params) for f in s.fruits)

        # Node-driven flowering: generic appearance gates set fraction
        if active && day >= gp.set_start_day && dnode > 0.0
            ss = dDM_pot / (veg_sink + fruit_sink_pre + 1e-9)
            set_frac = clamp((ss - gp.set_r_low) / (gp.set_r_high - gp.set_r_low), 0.0, 1.0)
            # Appearance rate scales new cohort per node increment (dnode drives timing)
            n_set = generic_appearance(Tmean_C, generic_params) * dens * dnode * set_frac
            n_set > 0.0 && push!(s.fruits, Fruit(0.0, 0.0, n_set))
        end

        # Allocate using M94 eq.4 / H96 eq.1 helper
        fruit_sinks_vec = [generic_sink(Tmean_C, f.dev, f.n, generic_params) for f in s.fruits]
        potentials = vcat([veg_sink], fruit_sinks_vec)
        growths, _ = allocate_daily(potentials, dDM_pot)
        veg_DM     = growths[1]
        new_leaf_area = veg_DM * gp.frac_leaf * gp.SLA
        new_leaf_area > 0.0 && push!(s.leaves, LeafCohort(new_leaf_area, 0.0))
        s.W_stem += veg_DM * gp.frac_stem
        s.W_root += veg_DM * gp.frac_root
        for (i, f) in enumerate(s.fruits)
            f.W += growths[i + 1]
        end
    end

    # --- harvest mature fruit --------------------------------------------------
    keep = Fruit[]
    for f in s.fruits
        if f.dev >= 1.0
            s.yield_DW += f.W
            s.yield_FW += (f.W / gp.DMC_fruit) / 1000.0
        else
            push!(keep, f)
        end
    end
    s.fruits = keep

    # --- natural senescence: oldest cohorts die past their lifespan ------------
    filter!(c -> c.age < gp.leaf_lifespan_dd, s.leaves)

    # --- management DE-LEAFING: trim the OLDEST leaf area to the working target -
    LAI = isempty(s.leaves) ? 0.0 : sum(c.area for c in s.leaves)
    if LAI > deleaf_target && !isempty(s.leaves)
        sort!(s.leaves, by = c -> c.age, rev = true)            # oldest first
        excess = LAI - deleaf_target
        i = 1
        while excess > 1e-9 && i <= length(s.leaves)
            if s.leaves[i].area <= excess
                excess -= s.leaves[i].area; s.leaves[i].area = 0.0
            else
                s.leaves[i].area -= excess; excess = 0.0
            end
            i += 1
        end
        filter!(c -> c.area > 1e-9, s.leaves)
    end
    s.LAI = isempty(s.leaves) ? 0.0 : sum(c.area for c in s.leaves)
    return nothing
end

"""
    simulate_cohort_measured(p, gp, season; LAI0, deleaf_target, generic_params) -> daily series

Cohort crop model driven by MEASURED greenhouse climate (like
simulate_growth_measured) for validation/calibration.

`generic_params`: optional `GenericParameters` instance passed through to
`grow_cohort!`. Default `nothing` uses the legacy V2.2 smoothstep sink.
"""
function simulate_cohort_measured(p::NamedTuple, gp, season; LAI0 = 0.5,
                                  deleaf_target = Inf, density_sched = nothing,
                                  layered = true,
                                  generic_params::Union{Nothing, GenericParameters} = nothing)
    # density_sched(day)::stems/m2 lets the grower change PLANT DENSITY in-season
    # ("separate the plants"): planting dense for fast canopy closure, then pulling
    # stems. A DROP on day d scales the standing leaf cohorts AND the unharvested
    # fruit by d_now/d_prev (the pulled stems take their leaves and green fruit with
    # them), and future leaf demand / fruit set use the reduced density. Default
    # nothing = constant gp.stem_density (identical to before).
    dsched = density_sched === nothing ? (_-> gp.stem_density) : density_sched
    d_prev = dsched(0)
    s = init_cohort_state(gp; LAI0 = LAI0 * d_prev / REF_DENSITY)
    day = Int[]; LAI = Float64[]; yieldFW = Float64[]; Wfruit = Float64[]
    Pg = Float64[]; node = Float64[]; nleaf = Int[]; dens_out = Float64[]
    for d in 0:(season.ndays - 1)
        d_now = dsched(d)
        if d > 0 && d_now < d_prev - 1e-9          # thinning event: pull stems
            r = d_now / d_prev
            for c in s.leaves; c.area *= r; end
            for f in s.fruits; f.n    *= r; end
            s.LAI *= r
        end
        d_prev = d_now
        M   = season.day_samples[d + 1]
        Pgd = layered ? daily_assimilation_layered(M, s.leaves, p, season.dt, gp.k_ext) :
                        daily_assimilation_measured(M, s.LAI, p, season.dt, gp.k_ext)
        grow_cohort!(s, Pgd, daily_mean_T(M), gp, d;
                     deleaf_target = deleaf_target, dens = d_now / REF_DENSITY,
                     generic_params = generic_params)
        push!(day, d); push!(LAI, s.LAI); push!(yieldFW, s.yield_FW)
        push!(Wfruit, isempty(s.fruits) ? 0.0 : sum(f.W for f in s.fruits))
        push!(Pg, Pgd); push!(node, s.node); push!(nleaf, length(s.leaves)); push!(dens_out, d_now)
    end
    return (; day, LAI, yield_FW = yieldFW, W_fruit = Wfruit, Pg, node,
              n_cohorts = nleaf, density = dens_out, state = s)
end
