# Cross-team growth calibration (AGC-2018, all 6 growers)

`test/calibrate_all_teams.jl` — full-length multi-objective t-walk (weekly yield +
leaf number + derived LAI) for every grower. Outputs `test/figures/teams_params.csv`,
per-parameter comparison plots `teams_<param>.png`, and per-team fits in `teams/`.

## Posterior means (full chains, 200k / 80k burn-in)

| team | J_max (µmol) | node_rate | veg_sink_max | LAI_max | set_start_day | yield model/obs |
|------|------|------|------|------|------|------|
| AiCU        | 119 | 0.092 | 24.4 | 2.50 | 16.3 | 32.8 / 37.3 |
| Croperators | 145 | 0.099 | 36.9 | 2.51 | 14.5 | 50.1 / 48.0 |
| DeepGreens  |  93 | 0.077 | 29.4 | 2.51 | 14.0 | 36.3 / 35.8 |
| Reference   | 178 | 0.093 | 35.9 | 2.51 | 16.4 | 48.6 / 49.4 |
| Sonoma      | 147 | 0.085 | 26.6 | 2.50 | 13.0 | 54.9 / 51.3 |
| iGrow       | 198 | 0.089 | 36.7 | 2.51 | 12.4 | 51.0 / 46.8 |

(J_max shown in µmol = value ×1e6; error bars tight, teams statistically distinct.)

## Assessment

- **Robust**: all 6 fit with the same model/priors; yields within ~10% (AiCU worst
  at −12%; most within a few %).
- **Development/timing CLUSTER** — `node_rate` 0.077–0.099, `set_start_day` 12–16 d.
  Consistent across growers ⇒ real, transferable physiology.
- **Source/allocation VARY ~2:1** — `J_max` 93–198 µmol, `veg_sink_max` 24–37;
  statistically real, all physically plausible, ordered ~by yield.

## Key caveat (confound)

The `J_max` / `veg_sink_max` spread is **partly confounded by the shared
plant-density assumption** (2.5 stems/m², 0.05 m²/leaf used for ALL teams). Since
those parameters scale with intercepted light per m² ground, a wrong per-team
density is absorbed into that team's `J_max`. `LAI_max ≈ 2.5` for all is **imposed**
by the same assumption, not independently identified. ⇒ absolute cross-team `J_max`
ranking is NOT clean physiology until real per-team densities are used.

## Recommendations

- Get per-team planting densities (AGC protocol docs, or a per-compartment value)
  and re-run — should sharpen `J_max`/`veg_sink` and likely tighten the spread.
- AiCU −12% is the outlier; its flush suppression may be slightly strong for its
  climate — a quick look if we want it tighter.
- Otherwise parameters are reasonable → ready for goal #2 (coupled ODE twin +
  RL control agent), keeping the density caveat in mind.

## Status vs project goal #1

Growth model (goal #1) is effectively complete: photosynthesis + climate physics
(V2.2, fast/type-stable) coupled to a Marcelis-style cucumber growth model,
calibrated and validated across 6 independent growers. Remaining polish: real
densities, optional harvest-flush amplitude tuning.
