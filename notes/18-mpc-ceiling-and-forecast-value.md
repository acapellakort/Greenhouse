# 18 — MPC perfect-foresight ceiling; the agent ignores its forecast

## The ceiling (test/mpc_ceiling.jl)
The twin is deterministic given a weather year, so the best OPEN-LOOP schedule
with the weather revealed is that year's optimum — an upper bound on any causal
controller. Approximated per TEST year with CEM (10-day blocks × 11 setpoints,
POP 48, 20 iters), warm-started from the deployed agent's realized trajectory.

    year   CEM-const   SAC   SAC-blockavg   MPC-ceil   gap SAC->ceil
    1990     -1.17     0.40      0.64          2.30        1.90
    1997      0.67     1.76      1.70          3.84        2.08
    2003     -1.50     1.20      0.76          2.74        1.54
    2009      0.33     1.99      2.19          3.88        1.89
    2014      1.94     2.98      3.65          5.41        2.43
    2015      0.44     1.60      1.98          4.00        2.39
    mean      0.12     1.66      1.82          3.70        2.04

Resolution check PASSED: SAC-blockavg (1.82) ≈ SAC (1.66) — 10-day blocks retain
the agent's behavior, so the ceiling is fair (not a coarsening artifact).

## forecast_skill semantics (IMPORTANT — earlier notes had this backwards)
In weather_episodes.jl `forecast_row`: `Tf = true + noise·sigma·skill`. So
**skill = 0 is a PERFECT forecast** (zero noise), skill = 1 nominal noise, larger
worse. The campaign eval env used `forecast_skill = 0.0` = PERFECT 3-day forecast
(NOT blind, as note 17 / earlier note 18 said). Training used skill = 1.0.

## Forecast-value result (test/forecast_value.jl) — the agent IGNORES the forecast
Deployed agent, mean TEST profit vs forecast noise:

    fcst_skill (0=perfect):  0.00   0.25   0.50   0.75   1.00
    mean TEST profit      :  1.66   1.63   1.60   1.63   1.63

FLAT. Even a PERFECT 3-day forecast (skill 0) gives no lift over a noisy one. The
agent's policy does not condition on the forecast at all. So SAC's entire
1.66-vs-0.12 advantage is REACTIVE (responding to current climate+crop state),
with ZERO anticipatory value. The ~2 EUR/m2 gap to the ceiling is timing/foresight
value the agent never learned to capture.

## Why the agent ignores the forecast (hypotheses)
1. Credit assignment: the payoff of an anticipatory action (pre-heat before a cold
   day) is delayed; 1-step TD may not propagate it against the strong immediate
   reactive signal.
2. Our warm-start fix (note 17) tethered the actor to the CEM CONSTANT program — a
   forecast-blind schedule — for the first 40k steps; this plausibly anchored the
   policy in a forecast-ignoring basin. The stability win may have cost forecast use.
3. Horizon: obs carries only a 3-day forecast (`_FC_H = 3`); the ceiling's value
   needs long-horizon timing, so even perfect 3-day info captures only a fraction.

## Caveat on the 3.70 ceiling
It assumes PERFECT 100-day foresight — physically impossible (weather is
unpredictable beyond ~1-2 weeks). The realistic operational ceiling is lower and
is best measured by RECEDING-HORIZON MPC with a realistic (degrading) forecast.
So 3.70 is a loose upper bound; the achievable target sits between 1.66 and 3.70.

## Fork for next step (decision pending)
A. Receding-horizon MPC controller (reuses the ceiling CEM engine): re-optimize
   each day over an N-day realistic forecast. Deployable, classical greenhouse
   control, and yields the REALISTIC operational ceiling. Highest tractable value.
B. Fix RL to use the forecast: longer horizon (_FC_H=7), soften/replace the
   constant-program tether with FORECAST-CONDITIONED demos (e.g. MPC-generated),
   so the agent learns anticipatory control. MPC-as-teacher synthesis.
C. Conclude the RL arc: reactive SAC (1.66) >> CEM (0.12); anticipatory headroom
   to the ceiling is future work. Write up CEM -> SAC -> forecast-value -> ceiling.

## Status
- Stomata concluded (note 16). inc1 + warm-start fix done (note 17): SAC 1.66 vs
  CEM 0.12, 14/15 seeds.
- MPC ceiling 3.70; agent ignores forecast (reactive-only). Direction pending.
