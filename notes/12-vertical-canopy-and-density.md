# 12 — Vertical canopy profile, cohort crop model, and the density question

Extends the cucumber growth model (V2.2 `GreenhouseSim`) with a **leaf-age cohort
canopy**, **node-driven thermal-age fruit cohorts**, and **layered (per-cohort)
photosynthesis** — the "vertical canopy profile" needed to reason about
de-leafing and plant density. Resolves what the earlier big-leaf de-leaf attempt
could not.

## What was added (files)

- `src/growth_cohort.jl`
  - `LeafCohort(area, age)` / `CohortState`: leaves are age cohorts produced at the
    head (node rate); senescence drops cohorts past `leaf_lifespan_dd`; de-leafing
    trims the **oldest** (bottom, shaded) cohorts to a working target.
  - Node-driven flowering: each node's flowers set (carbohydrate-regulated) into
    `Fruit` cohorts carrying `dev` = thermal age; harvested at `dev ≥ 1`.
  - **`daily_assimilation_layered`** + `canopy_light_profile`: integrate the FvCB
    leaf rate **down the light gradient** `A(I₀·e^{-k·L})` slab by slab
    (age-sorted, midpoint rule), instead of running it at top light × `lai_eff`.
  - `simulate_cohort_measured(...; layered=true, deleaf_target, density_sched)`:
    `density_sched(day)→stems/m²` models "separate the plants" — a mid-season
    density drop scales standing leaf **and unharvested fruit** by the ratio.
- `configfiles/cucumber_growth.json`: calibrated cohort leaf production
  `leaf_per_node=0.040`, `leaf_lifespan_dd=600` → natural LAImax ≈ 2.5 (matches
  identified `LAI_max`). Added `LAI_deleaf`.
- Tests: `test/test_cohort.jl` (leaf-production sweep + de-leaf),
  `test/test_density.jl` (big-leaf vs layered, LAI optimum, de-leaf, schedule),
  `test/test_density_econ.jl` (k_ext × density, profit proxy).

## Key results (AiCU measured climate, season 2018)

1. **Layered photosynthesis fixes *more* carbon than big-leaf — textbook-correct.**
   Big-leaf ran the leaf rate at top-of-canopy (light-saturated, low per-photon)
   light and scaled by interception; the layered integral lets shaded leaves work
   at higher quantum efficiency. Total Pg rises: 2.5 stems 26.1→33.5 kg,
   4.0 stems 29.0→39.1 kg.
2. **The ~15% cohort yield gap was mostly the big-leaf light artifact, not
   turnover.** At natural density (2.5 stems) the layered cohort model gives
   **33.5 kg ≈ measured ~31** with no source re-calibration.
3. **Optimal density is governed by `k_ext` (light extinction).** Yield is
   single-peaked; the peak walks down as extinction rises:

   | k_ext | optimal density | yield (kg/m²) |
   |------:|:---------------:|:-------------:|
   | 0.60  | 5.0 stems | 48.9 |
   | 0.75  | 4.5       | 39.3 |
   | 0.90  | 3.5       | 32.7 |
   | 1.10  | 3.0       | 26.7 |

   `k_ext=0.90` also reproduces measured yield at a plausible 2.5–3 stems, hinting
   the true value is ≈0.85–0.95, not the assumed 0.75.
4. **Economics barely shift the optimum.** Plant + CO₂ marginal cost ≈ 2.5 €/m² vs
   ≈ 35 €/m² revenue, so profit-optimal density ≈ yield-optimal (4.5 at k=0.75,
   3.5 at k=0.90). Density is an agronomic/light decision, not a cost decision.
5. **Dense planting to the light-limited optimum is confirmed** (denser → higher
   yield up to the peak), **but separating the plants later never pays** in the
   carbon+profit balance: thinning discards standing fruit, and with a correct
   light profile the retained deep leaves are not dead weight. De-leaf (leaf only)
   is at best neutral.

## The confound (important)

Plant density is **not recorded** in the AGC-2018 AiCU dataset (fixed protocol
value; derived LAI rides on assumed `leaf_area`, `N_window`). `k_ext` and
`stem_density` both scale interception, so **yield alone cannot separate them** —
many (k, density) pairs reproduce ~31 kg. Density is therefore a **control
variable the agent sets**, and `k_ext` should be fixed from canopy architecture
(literature) rather than fit against this data.

## Implications for the control agent (goal #2)

- Density lever = **choose the optimal constant density**, set by `k_ext`
  (≈3–3.5 stems at k≈0.9). A dense-early/thin-late *schedule* is not
  carbon/profit-justified in the current model.
- The AGC challenge simulator *does* express density as a declining schedule
  (e.g. `@plantDensity "1 56; 25 45; 50 35; 65 20"`, from a v2.0-challenge scratch
  run) — so growers/organizers do thin. That benefit is **non-carbon**
  (humidity/Botrytis, labour, fruit grade) and would require adding a
  humidity/disease penalty to appear in the model.

## Open items

- Fix `k_ext` for AiCU high-wire cucumber from literature/leaf-angle (≈0.85–0.95).
- (Optional) add a humidity/disease penalty so dense-canopy downsides are
  represented and thinning can pay for the right reasons.
- Integrate the cohort crop model into the V3.0 twin (`simulate_twin` still uses
  big-leaf `grow!`) so net profit under density scenarios is metered end-to-end.
- Then resume the RL agent build (Stage 3): the density action can be a single
  season-constant choice rather than a per-step schedule.
