# 17 — inc1 campaign, training-instability diagnosis, warm-start fix (RESOLVED)

## The inc1 actuator increment
Added to the RL control surface (twin unchanged, constant-g): continuous **lamp
dimming** (U12 as a 0–1 intensity, the 11th action) and **flue-gas CO2** (boiler
combustion supplies free CO2 while heating: `paid = max(dose − 0.85·eta13·heatW, 0)`).
Density set to the inc1 CEM optimum 2.40. `WARM_PHYS` = inc1 CEM program.
Env/CEM/warm-start all at DENSITY=2.40 for a fair comparison.

## inc1a: broken (60k steps, no warm-start protection)

    CEM constant  : TEST 0.12
    SAC best-ckpt : TEST -0.13 +/- 0.40 ; beats CEM 4/15 seeds
    deployed (seed 112): TEST 0.59

Robustness collapsed vs the earlier 10-action campaign (+1.48±0.68, 15/15).

## Diagnosis: truncated recovery + destructive early exploration
Per-seed VAL curves all had the same shape: warm-start ≈ CEM → training DESTROYS
it (VAL crashes to −1.5..−2.5 by ~18–27k as the actor bolts from CEM against an
untrained critic; alpha collapses to the 0.010 floor) → slow recovery crossing
CEM only near 50–60k, if it finishes. Several seeds still climbing at the 60k
cutoff → recovery truncated. Root cause: actor updates start against a garbage
early critic. The 11th action + flue-gas coupling pushed the phase transition to
~50k, and 60k truncates. Evidence-based answer to "bigger campaign / longer runs
/ architecture": longer runs help; architecture does not (good seeds reach good
policies); more seeds secondary.

## Fix applied (train_sac.jl + train_sac_long.jl)
1. **Critic-only warmup**: `sac_update!` gained `update_actor::Bool`; train_one
   runs `CRITIC_PRETRAIN = 4000` critic-only updates on the demo-anchored buffer
   before the actor moves, so Q is accurate around CEM first.
2. **Decaying BC tether**: `sac_update!` gained `bc_target`, `bc_weight`; actor
   loss adds `bc_weight · mean((μ − atanh(2·WARM_A−1))²)`. `bc_weight` decays
   `LAMBDA_BC = 3.0 → 0` over `BC_DECAY = 40000` steps.
3. **Longer runs**: `total_steps 60k → 100k`, `eval_every 4000`, `BUDGET 8h`.

Smoke check confirmed the mechanism: step-500 VAL went from −32.6 (old) to +0.14
(new) — the early bolt is gone.

## inc1b: FIXED (100k steps, warm-start protection) — 15 seeds, 6.4 h

    CEM constant  : VAL 0.12  TEST 0.12
    SAC best-ckpt : TEST 0.76 +/- 0.58 ; beats CEM 14/15 seeds ; best 1.71, worst -0.11
    per-seed TEST : 0.36 1.71 0.38 1.66 0.22 0.35 0.98 0.92 1.12 1.69 0.46 -0.11 0.37 0.91 0.45
    deployed (best-VAL seed 104): TEST 1.66 vs CEM 0.12 -- beats CEM on ALL 6 years:
      year   CEM    SAC        year   CEM    SAC
      1990  -1.17   0.40       2009   0.33   1.99
      1997   0.67   1.76       2014   1.94   2.98
      2003  -1.50   1.20       2015   0.44   1.60   (mean 0.12 -> 1.66)

The fix worked: 4/15 -> 14/15, mean TEST −0.13 -> 0.76. The **deployed agent
(seed 104) beats CEM on every held-out year**, advantage +1.54 EUR/m2 — matching/
exceeding the old 10-action +1.48. The richer actuator set (lamp dimming +
flue-gas CO2) is now productively exploited. Best actor: test/sac_actor_best.jls.

Caveats / open: mean-over-seeds advantage (+0.64) is below the old 10-action mean
(+1.48), but environments differ (density 2.40 vs 2.90, different CEM baseline,
11 vs 10 actions) so it is not apples-to-apples. One seed (112) still lost
(−0.11) — likely slight undertraining; not investigated. Whether VAL had fully
plateaued at 100k not checked.

## Status
- Stomata increment concluded (note 16: autumn light-limited, flag OFF).
- inc1 (lamp dimming + flue-gas CO2) + warm-start fix: DONE, strong result.
- NEXT: MPC perfect-foresight ceiling — gives the SAC(1.66)-vs-CEM(0.12) gap a
  denominator (how close to the achievable optimum is the agent?).
- Optional: diagnose_policy on seed 104 to confirm it exploits the new levers
  (lamp dimming timing, cheap-CO2-while-heating).
