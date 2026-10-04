# 14 — Cohort crop in the twin, lamp-blindness fix, and RL vs CEM

Integrated the layered cohort crop (note 12) into the RL environment and re-ran
both controllers. Two substantive twin fixes and — after a multi-seed check — a
sober RL-vs-CEM conclusion.

## What changed (`rl_env.jl`)

`TwinEnv` gained a `crop::Symbol` switch (`:bigleaf` default, `:cohort`):
- `:cohort` uses `CohortState`/`grow_cohort!` with **layered per-cohort Pg**
  computed from the day's SIMULATED climate (`_pg_matrix` reconstructs hourly
  [Tair; CO2; canopy light] from the ODE solution, then `daily_assimilation_layered`).
- Handles episode-level density and mid-season thinning.
Both CEM (`optimize_cem.jl`) and SAC (`train_sac.jl`) reach the crop through the
env, so `crop=:cohort` switches the whole agent stack. Big-leaf path retained.
Action space extended to **10** (added `light_start` = lamp placement). An
**ODE-instability guard** was added to `env_step!`: an incomplete/non-finite solve
returns a −3 €/m² penalty and preserves the last good state.

## Twin fix: the big-leaf twin was lamp-blind

Big-leaf photosynthesis used `I9 = (1-eta1)·tau1·eta2·Idocel` — **solar only, no lamp
term** — so lamps cost electricity+heat but grew nothing, and every optimizer turned
them off (big-leaf 18 h-lamp: 13.7 kg ≈ unlit 13.2, −54 €/m²). The cohort `_pg_matrix`
**includes lamp PAR** (`alpha12·Light_on`), fixing this.

## The environment came alive (cohort crop)

Cohort crop, same programs (1995/2007): COLD (no lamps) 11 kg ≈ −1.7 €/m²; LIT (18 h)
**33 kg**, −35 €/m². Neither extreme optimal ⇒ a genuine **interior lamp optimum**.
**CEM on the cohort twin**: −24.1 → **−0.77 €/m²**; light_hours 3.19, density 2.90,
yields 18–20 kg; the ~3 h from 04:00 sit mostly **off-peak** (<07:00).

## RESULT (corrected by multi-seed): a tuned constant program is the robust controller

A **single** warm-started SAC seed (pre-guard code) looked like a clear win — held-out
−0.74 vs CEM −2.19, worst-year −35.9 vs −46.2. **But a 3-seed replication (guarded
code, `train_sac_multiseed.jl`) does NOT reproduce it:**

| controller (cohort twin) | tuned 4 yr | held-out 27 yr | all 31 | worst yr |
|--------------------------|-----------:|---------------:|-------:|---------:|
| CEM (constant)                        | −0.80 | **−2.19** | −2.01 | **−46.2** |
| SAC single favorable seed (pre-guard) | −0.23 | −0.74 | −0.68 | −35.9 |
| **SAC 3-seed mean ± std (guarded)**   | −0.63 ± 0.89 | **−2.56 ± 0.21** | −2.31 ± 0.29 | **−55.5 ± 10.3** |

Per-seed held-out: −2.79, −2.36, −2.55. **SAC beats CEM on held-out in 0/3 seeds**,
is worse on average (−2.56 vs −2.19), and is *less* robust on the worst year
(−55 vs −46). The seeds cluster tightly at ≈−2.5 and *degraded* from the warm-start's
≈−1.2 start — i.e. SAC training on this near-static-optimum problem is high-variance
and, on average, does not improve on the best constant program.

**Conclusion:** on this twin, direct optimization of a compact constant/MPC-style
program (CEM) is the robust, pragmatic controller; model-free deep RL (SAC) is
high-variance and does not reliably beat it. Consistent with the greenhouse-control
literature (MPC / direct optimization competitive with or superior to model-free RL;
AGC winners blend RL+MPC). The earlier single-seed "win" was a favorable draw and,
being pre-guard, not even the same code path.

## Mechanism of the best seed (not a reliable advantage)

`diagnose_policy.jl` on the favorable seed showed **phenological light-budget
allocation** (not weather-timing: `corr(lamp_hours, solar) ≈ 0`): lamps ≈0 during
establishment (days 0–25) and wind-down (80–100), 3–8 h in the productive window
(30–70), at a *lower* mean (1.7–2.3 h) than CEM's 3.19 h. This is a sensible,
genuinely state-dependent strategy — but the multi-seed result shows SAC does not
learn it *reliably*, so it is not a dependable advantage over the constant program.

## Follow-ups
- If pursuing RL further: reduce variance via best-checkpoint selection on a held-out
  *validation* split, lower LR, larger networks, or more seeds; and re-run the
  favorable config under the guard to isolate the guard's effect.
- Pragmatic path: adopt the CEM/MPC constant program; optionally a receding-horizon MPC
  over the twin as the deployable controller.
- Model documentation: `docs/greenhouse_model.{tex,pdf}`; agent methods/results:
  `docs/agent_results.{tex,pdf}` (being corrected to the multi-seed conclusion).
