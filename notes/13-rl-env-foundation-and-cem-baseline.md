# 13 — RL environment foundation + CEM baseline + SAC (Stage 3b)

Built and validated the pieces the RL agent needs, a gradient-free baseline
("number to beat"), and a working SAC agent. All in V3.0 `GreenhouseControl`.

## VPD humidity control (was RH)

The humidity actuator now triggers on a **VPD floor** (kPa), not relative
humidity. `rhs_control!` computes `VPD = (Pws(T2) − V1)/1000` and dehumidifies
(crack energy screen to 0.85, add roof/side venting proportionally) when
`VPD < VPDmin` (default 0.5 kPa). Setpoint field `RHmax → VPDmin` through the
struct, `daynight_setpoints`, `screen_control`, `humidity_vent`, `heat_audit.jl`.
VPD is the grower-standard variable (AGC used `HumDef`); gives the agent a
temperature-independent humidity signal. Backup: `climate_computer.jl.bak_prevpd`.

## Weather episodes + noisy forecast (`weather_episodes.jl`)

- `WeatherBank(csv; plant_month, plant_day, ndays)` parses the 32-year Holland
  file **once** and lists fully-covered start-years. Result: **31 usable years**
  (1988–2017, 2023); the gap years **2018–2022 are correctly excluded**.
  `weather.jl` refactored into `prepare_meteo` + `weather_from_df` (one source of
  truth; loud error on empty window) so episodes slice the cached frame, not re-read.
- `SeasonForecast` / `forecast_row(day)`: noisy H=3-day daily forecast of the
  control-moving drivers (outside T, radiation, wind). Lead-dependent error
  (T ±1→2.5 °C, rad ±15→35%, wind ±1→2 m/s) × a **skill knob** (0 = perfect …
  large = useless) for the forecast-value ablation. Validated: skill 0 reproduces
  truth; lead-3 spread grows with skill.

## Stepping env (`rl_env.jl`)

- **Action = 9 daily setpoints** (normalized): Tday, Tnight, CO2, light_hours,
  VentpBand, ofset, ToutMax, **VPDmin, ScreenRad**. Plant **density** is an
  episode-level `density_sched` set at reset. Reward = daily net profit €/m².
- **Observation = 18**: 6 crop (day, LAI, yield, node, leaf, fruit) + 3 climate
  summary (yesterday T, CO2, VPD) + 9 forecast (3 days × Tout/rad/wind).
- `TwinEnv(params, gp, bank; ...)` samples a fresh covered year each reset
  (domain randomization); fixed-weather constructor retained.

## CEM baseline (`optimize_cem.jl`, threaded)

Cross-Entropy Method over a **constant-season program (9 setpoints + density)**,
scored as mean net profit over 4 fixed years (1995/2001/2007/2013).

**Result: −43.79 → −0.29 €/m² (Δ +43.5), ~break-even, monotone over 12 iters.**
Winning program: **lamps 0 h**, **T 16 °C day+night (lower bound)**, CO2 790,
VPDmin 0.20, density **2.29 stems**. Per-year −0.21 / +0.23 / −1.10 / −0.07.

### Reward audit (important)
- The optimum is **cost-dominated** and sits on the **temperature lower bound**:
  22→16 °C + no lamps costs ~1 kg yield but saves ~€43. The big-leaf twin's yield
  is light/source-limited over the season, so temperature barely moves total kg —
  "grow cold and cheap" is the true optimum of this reward.
- **Tight per-year spread (~1.3 €/m²)** ⇒ a reactive (forecast-driven) policy's
  headroom over a constant program is **modest on this twin**. The big lever is
  swapping in the **cohort crop** (note 12; ~31 kg) so revenue dominates and
  warmth/CO2 become worth paying — then the control tradeoffs get rich.

## SAC agent (`test/train_sac.jl`, Julia-native Flux)

Hand-rolled SAC (kept OUT of the core package so non-RL runners don't pull Flux):
twin Q-critics + squashed-Gaussian actor (tanh log-prob correction) + auto entropy
temperature + replay buffer + soft target updates via `Functors.fmap`. Reward scaled
×10 for learning (reported euros stay true); density fixed at 2.3 (CEM optimum) to
isolate the climate policy. Requires `Pkg.add("Flux")`.

**Training (24k steps, UTD 2, domain-randomized weather, noisy forecast skill=1):**
monotone, stable convergence of deterministic-eval profit:
−6.74 (4k) → −3.58 (8k) → −2.77 (12k) → −2.45 (20k) → **−2.13 €/m² (24k)**;
α collapses 0.18 → 0.008 (policy near-deterministic).

**Read:** the full RL pipeline works end-to-end. SAC lands €1.8/m² behind CEM's
−0.29 — expected, not a defect: CEM overfits a fixed program to the 4 eval years,
while SAC learns a state-dependent policy across all 31 years under a *noisy*
forecast (a harder, more general problem). Both near break-even because the twin is
cost-dominated (reward audit). SAC's advantage should show as **generalization /
robustness**, and its **ceiling rises with the cohort crop**.

## Next (options)
- **Generalization test**: score the CEM fixed program vs the SAC policy on ALL 31
  years (and/or held-out years) — the real justification for a learned policy.
- **Integrate the cohort crop** (note 12) into `simulate_twin` so revenue dominates
  and the control tradeoffs (warmth, CO2, density, de-leaf) become economically live.
- Hand **density / de-leaf** to the agent as decisions; warm-start SAC from CEM to
  close the small gap; forecast-skill ablation (is the forecast worth anything?).
