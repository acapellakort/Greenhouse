# 10 — Stage 3: RL agent design, economics reward, and weather/forecast

Goal #2: an AI agent that sets the greenhouse climate program to maximize net
profit, trained on the validated V2.2/V3.0 twin. This note records the design
decisions and the Stage-3a foundation (resource accounting + reward), built and
validated.

## The environment (MDP)

- **State/observation**: crop (LAI, biomass, node, fruit load, day-of-season),
  climate summary (recent air T, CO2, RH), and a WEATHER FORECAST (below).
- **Action**: setpoints. Start with 3 core levers (day temp, night temp, CO2)
  chosen daily; optionally the full program (vent band, offset, light hours, ToutMax).
- **Reward**: daily net profit €/m² (below).
- **Transition**: one day of `simulate_twin` (~2 s/season — fast, deterministic).
- **Episode**: one ~100-day season; a different weather year sampled per episode.

## Algorithm landscape (state of the art for this problem)

Because the twin is fast, deterministic and differentiable, more than model-free
RL applies:
- **Model-free continuous control**: SAC (sample-efficient, robust — first choice),
  TD3/DDPG, PPO (robust baseline, sample-hungry but the sim is cheap). Small MLPs
  suffice; add a temporal encoder only if feeding a forecast sequence.
- **Exploits our model**: MPC / receding-horizon over the twin (strong, interpretable
  baseline); differentiable-simulator policy gradient (backprop season profit through
  the ODE via SciMLSensitivity — Julia-native); gradient-free (CMA-ES/CEM) over a
  parameterized policy (simplest first baseline).
- Greenhouse-RL literature centers on economic-reward DRL (PPO/SAC), reward-component
  auditing (economic rewards are easy to game), and grower-in-the-loop RL. AGC winners
  mix RL and MPC.

**Plan**: (1) accounting + reward + non-learned baseline (MPC/CMA-ES) to get a
number to beat and to audit the reward; (2) SAC as the learned agent, trained across
weather years; (3) PPO cross-check. Ecosystem: keep the fast Julia twin; either Julia
RL (ReinforcementLearning.jl/Crux.jl) or wrap it Gym-style for Python SB3/CleanRL
(more mature) — user has no preference, decide when building the env API.

## Net-profit reward (economics.jl, built)

Coefficients from note 08; cucumber price 0.889 €/kg (AGC-2018 Reference
Prod_value_cum / Total_Prod_cum). Net = revenue − heating(0.09 €/kWh) −
CO2(0.30 €/kg) − electricity(0.30 peak / 0.20 off-peak €/kWh) − fixed
(greenhouse 15 €/m²/yr prorated + plant cost). `simulate_twin` now meters
per-day heat_kWh, co2_kg, lamp_kWh_peak/off (`day_resources`, mirrors the
controller actuator formulas); `season_economics` prices them.

**Validation vs the real AGC-2018 Reference grower (per day):** CO2 0.088 vs 0.079
kg/m²/d ✓; lamp elec 1.18 vs 1.26 kWh/m²/d ✓ — the accounting is trustworthy.

## Realism audit and the FIR/sky fix (done)

Heating first came out ~1.8× high (2.12 vs 1.18 kWh/m²/d). A season **heat-balance
audit** (V3.0 `test/heat_audit.jl` — re-derives every dT2 term and integrates it;
consistency-checked to ~1% once sampled at 5-min) found the cause: **FIR-to-sky was
the dominant loss, 394 kWh/m² (~164 W/m²)**, because the twin radiated to a fixed
clear −0.4 °C sky every night. NOT the floor (floor on/off heating 212 vs 214) and
only weakly the night setpoint (19→17 °C: −11%).

Fix: `load_weather(...; sky_from_clouds=true)` computes the effective sky temperature
per hour from air T, RH and cloud cover (Berdahl-Martin clear-sky emissivity → Tout as
cloud rises). Result: FIR 394→230, **heating 208→129 kWh/m² (1.29 vs the reference
1.18)** — over-heating resolved. All V3.0 runners now pass `sky_from_clouds=true`.
Documented in FROZEN.md; consistent with the venting calibration (measured sky).

Remaining (not blocking): the twin's absolute **yield** (19 kg/m², 1998 summer) is
unverified against a matched season — the 32-year file lacks 2018 Aug-Dec, so a fair
twin-vs-49 kg check needs the AGC-2018 weather in a load_weather-compatible form. The
growth model itself was validated to ~49 kg on measured climate, so this is a
season/weather comparison, not a known defect.

## Weather & forecast (design)

- **Training weather**: `dataset_meteo_holanda.csv` — hourly Holland, nominally
  1988–2023 but WITH GAPS (2018 ends mid-July, 2019 missing). The RL env must FILTER to
  years with complete season coverage before sampling one per episode (domain
  randomization). `load_weather` should also fail loudly on an empty date range.
- **Forecast given to the agent**: only the drivers that move control — outside
  temperature, solar radiation, wind (not humidity/CO2). Horizon 3 days, aggregated to
  daily values; beyond that, climatology.
- **Uncertainty**: forecast = truth + per-issuance correlated error growing with lead
  time — temp ±1 °C(d1)→±2.5 °C(d3), radiation ±15%→±35%, wind ±1→±2 m/s. Add a skill
  knob (0 = perfect … large = climatology) for the "how much is the forecast worth"
  ablation. Regenerate the noisy forecast each episode; the twin runs on the true
  sampled-year weather.
