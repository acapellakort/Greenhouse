# Cucumber growth model — design & coupling (Marcelis source–sink)

Goal: add a plant/vegetable growth model to the V2.2 digital twin, closing the
loop photosynthesis → dry matter → LAI → back into climate + photosynthesis.
Crop is **cucumber** (dataset = Autonomous Greenhouse Challenge 2018; Production.csv
gives cumulative cucumber fresh yield, ~36 kg FW/m² over the season).

Chosen approach: **Marcelis-style source–sink dry-matter partitioning** (Marcelis
1994; Marcelis et al. 1998), driven by the existing FvCB assimilation `A`, on a
**daily** timestep — the outer loop the 2.0 README always intended for growth.

## Time structure (two nested loops)

- **Within a day (seconds):** the fast climate+FvCB ODE runs at a *fixed* LAI,
  and assimilation `A` is integrated over the 24 h → daily gross assimilate `Pg`.
- **Once per day:** the growth model spends `Pg` on maintenance + growth, sets and
  grows fruit cohorts, harvests mature fruit, prunes leaves, and updates biomass
  and LAI. The new LAI feeds the next day's within-day ODE.

This also fixes the 2.0 `LAI = I9` quirk: the growth path uses the real LAI for
both canopy energy balance and photosynthesis.

## State variables (per m² ground, dry weight g)

- `W_leaf`, `W_stem`, `W_root` — vegetative pools
- `LAI = SLA · W_leaf`
- `fruits` — cohorts, each `(dev∈[0,1], W [g/m²], n [fruits/m²])`
- `yield_DW`, `yield_FW` — cumulative harvested fruit
- `node` — main-stem node number (development clock for fruit set)

## Daily update (Marcelis source–sink)

1. **Gross assimilate** `Pg` [g CH₂O m⁻² d⁻¹] = (∫A dt over day, mg CO₂ m⁻²)
   ·1e-3·(30/44).  (A integrated as an extra ODE state.)
2. **Maintenance respiration** `Rm = Q10^((T̄−25)/10) · Σ k_m,i W_i`.
3. **Assimilate for growth** `Cav = max(Pg − Rm, 0)` → `ΔDM_pot = Cav / asrq`.
4. **Development**: `node += r_node(T̄)`, each cohort `dev += 1/D_fruit(T̄)`.
5. **Fruit set**: add a cohort with `n = set_rate(T̄, node)` fruits.
6. **Sink strengths** (potential growth rates, g DM d⁻¹):
   - per fruit: `Wf_max · 6·dev·(1−dev) / D_fruit_days` (smoothstep integral = Wf_max)
   - vegetative: `veg_sink(T̄) · (1 − LAI/LAI_max)` (vegetative demand falls as canopy fills)
7. **Partition** available DM by sink share: `ΔW_i = (sink_i / Σsink) · ΔDM_pot`
   (source-limited). Vegetative split leaf/stem/root by fixed fractions.
8. **Harvest**: cohorts with `dev ≥ 1` → `yield_DW += W`, `yield_FW += W / DMC_fruit`, drop.
9. **Leaf pruning**: hold `LAI ≤ LAI_max` (remove oldest leaf DW to senescence).
10. `LAI = SLA · W_leaf`.

Performance note: the daily step runs ~100×/season, so it uses ordinary Julia
(vectors of cohorts, a mutable struct) — only the within-day ODE must be fast, and
that path is unchanged and already type-stable.

## Parameters (cucumber literature defaults — all calibratable)

`SLA≈0.030 m²/g`, `k_m_leaf≈0.03, k_m_stem≈0.015, k_m_root≈0.01, k_m_fruit≈0.01
g CH₂O g⁻¹ d⁻¹`, `Q10_resp=2`, `asrq≈1.4 g CH₂O/g DM`, `frac_root≈0.08`,
`leaf:stem split ≈ 0.65:0.35` of vegetative, `LAI_max≈2.5–3`, `T_base≈10 °C`,
`D_fruit≈12 d` (fruit growth period), `Wf_max≈16 g DW/fruit` (~400 g FW · 4% DM),
`DMC_fruit≈0.04`, `set_rate` ~1 fruit m⁻² d⁻¹ (tunable), `veg_sink` tunable.
These go in `configfiles/cucumber_growth.json`, loaded into a NamedTuple exactly
like the climate/crop params, so any of them can be inferred later.

## Calibration data (from the AGC-2018 dataset, per team folder)

- `Production.csv`: `Total_Prod_cum` (kg FW/m²), `ProdA/B_num` (fruit count) → the
  primary yield target (compare model `yield_FW`).
- `CropManagement.csv`: `N_leaves`, `LeafFormRate`, `Stem_elong`, `FruitGrw`,
  `Pruning` → constrain node/leaf dynamics and pruning.
- `Greenhouse_climate.csv`: `Tair`, `CO2air`, `RHair` → drive/validate the climate side.

## Coupling implementation (V2.2)

- `SimContext` gains an `lai` field (5-arg constructor defaults it to `NaN`, so the
  validated climate `rhs!` and `forward_map` are untouched).
- New `rhs_dyn!`: 6-state ODE `[T1,T2,C1,V1,ie,∫A]` that uses `ctx.lai` for the
  canopy (I1) and for photosynthesis (real LAI, quirk fixed), accumulating `A`.
- `simulate_growth(...)`: the season driver running the two nested loops.
- Cucumber params in `configfiles/cucumber_growth.json`; demo in `test/run_growth.jl`.

## Status / caveats

- This is a **v0**: structurally complete and calibration-ready, but the parameter
  values are literature placeholders — expect to fit `set_rate`, `veg_sink`,
  `Wf_max`, `SLA`, `asrq` against `Production.csv` before trusting yields.
- Fruit set / node dynamics are simplified vs full Marcelis (no explicit
  node-by-node cohort of leaves); can be refined if the yield/LAI fit needs it.
