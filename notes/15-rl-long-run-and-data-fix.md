# 15 — RL long-run campaign: SAC reproducibly beats CEM; demo-anchored replay lifts it

## Headline (best result): demonstration-anchored replay

On clean 30-year data, disjoint TRAIN/VAL/TEST splits, validation best-checkpointing, and
**demonstration-anchored replay** (25% of each minibatch drawn from a protected CEM/expert
buffer), SAC beats the tuned constant (CEM) program on held-out weather in **every seed**:

| controller (cohort twin, TEST = 6 held-out yr) | mean TEST €/m² | std | beats CEM | worst seed |
|-----------------------------------------------|---------------:|----:|----------:|-----------:|
| CEM constant                                  | −0.37 | — | — | — |
| SAC, plain replay (15 seeds)                  | +0.40 | 1.12 | 12/15 | −1.02 |
| **SAC, demo-anchored replay (15 seeds)**      | **+1.48** | **0.68** | **15/15** | **+0.55** |

Demo-anchored replay ~4x the mean (+0.40→+1.48), cut the std (1.12→0.68), took every seed
above CEM (12/15→15/15), and lifted the worst seed from −1.02 to **+0.55** (all seeds now
profitable). Best seed 115: TEST +2.89, beats CEM on all 6 TEST years (+1.08…+4.50).
Per-seed TEST (demo): 2.22,1.49,0.55,1.69,0.63,1.34,0.91,1.17,1.93,2.13,1.50,0.90,2.10,0.70,2.89.
Mechanism: the protected expert buffer keeps the critics grounded, so the late "phase
transition" fires reliably and lands on a good plateau in every seed rather than only ~12/15.

## The reproducible SAC-beats-CEM result (plain replay, for the record)

15 seeds, best-val ckpt: TEST **+0.40 ± 1.12**, 12/15 beat CEM, best +3.26. This already
REVERSED notes 13–14 ("SAC can't beat CEM"), which were flawed by (a) too-short training
stopping before the phase transition, (b) selection on the tuning years, (c) the corrupt
2023 year in the pool.

## Why the earlier runs failed: phase transition + rigorous evaluation

Every seed's validation curve: early exploration **crash** (−5…−33), long plateau (~−1.5),
then a **phase transition ~step 30–33k** climbing to +2…+3. Earlier 24k-step runs stopped
before it. The four fixes that flipped the verdict: clean data; held-out TEST never used for
selection; validation best-checkpointing; enough steps + gentle exploration
(lr 1e-4, target-entropy −3, no random warmup) + stiff-solver fallback in `env_step!`.
**The methodology was the finding, not the algorithm.**

## Data-quality bug fixed (2023 all-zero year)

`dataset_meteo_holanda.csv` year 2023 has full row count but all-zero temperature/radiation
(placeholder) → ~−47 €/m² for every policy; it had polluted the "worst-year" numbers in
notes 13–14 (disregard those). `valid_start_years` now rejects degenerate windows
(`max(rad)>50`, `std(T)>0.5`); usable years 31 → 30.

## Mechanism CONFIRMED on held-out TEST years (`diagnose_policy.jl`)

Not weather-timing (`corr(lamp_hours, solar) ≈ 0`). **Phenological light-budget allocation**:
lamps ≈0 in establishment (days 0–20), 5–7 h through the productive window (30–70), ≈0 in
wind-down (80–100); lamp **start placed early-morning off-peak** (<07:00 tariff); **CO2 ramped
up (400→1100+) with the maturing canopy**; temperature at the 16 C cost floor except an
establishment warm-push. All 4 checked TEST years profitable (+2.03/+2.66/+4.79/+3.40).

## Files
- `test/train_sac_long.jl` (splits, best-ckpt, demo-anchored replay via `DEMO_FRAC`,
  per-seed actor saving + `campaign_manifest.csv`); `test/sac_actor_best.jls`,
  `test/campaign_actors/`. `test/ensemble_eval.jl` (top-k action-average → TEST).
- Data fix in `src/weather_episodes.jl`; stiff fallback in `src/rl_env.jl`.

## Done / next
1. (DONE) Mechanism confirmed.
2. (DONE) Demonstration-anchored replay — big win (above).
3. **Ensemble** — `ensemble_eval.jl` top-5, pending (expect ≥ best single, better worst year).
- Setup context: autumn crop only (planted Aug 14, the hardest/lowest-light slot — 3 NL
  planting seasons exist); economics sells FRESH weight (yield_FW = dry/DMC_fruit, DMC=0.04)
  at 0.889 €/kg FW; tariffs are AGC-2018/2024 (pre-energy-crisis, likely low) — pending update.
- Actuator-freedom candidates (ranked for autumn): 2nd screen, continuous lamp dimming,
  min-pipe-temp, flue-gas CO2, richer temp schedule, CHP+buffer+electricity selling.
- `docs/agent_results.{tex,pdf}` to be updated with the demo-anchored + ensemble numbers.
