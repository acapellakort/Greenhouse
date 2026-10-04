# 19 — MPC -> RL distillation + confirmations (forecast, AGC, price, OOD)

## Idea
Reactive SAC ignored its forecast (note 18): 1.66 TEST, ~45% of the perfect-
foresight ceiling (3.60-3.70). Fix: teach anticipation by IMITATING a foresightful
expert. Generate forecast-driven expert demos with the open-loop CEM optimizer
(perfect foresight per TRAIN year, warm-started from SAC -> gen_mpc_demos.jl, 1800
transitions, mean expert profit ~2.3), then retrain SAC anchored to them with
STATE-CONDITIONED behavior cloning instead of the forecast-blind constant tether
(train_sac_mpc.jl). Demo = (obs w/ realistic noisy forecast) -> (action chosen with
full foresight). Privileged-teacher distillation.

## Result (100k steps, 15 seeds, 6.7 h)
    CEM constant : TEST 0.12
    reactive SAC : TEST 1.66
    SAC-MPC      : TEST 2.35 +/- 0.30 ; beats CEM 15/15 ; best 2.83 worst 1.73
    deployed (seed 202): TEST 2.83, beats CEM ALL 6 years, mean 0.12->2.83
    perfect-foresight ceiling: ~3.65
Deployed 2.83 = +70% over reactive 1.66, ~78% of the ceiling (was 45%), std 0.30.
Best actor: test/sac_actor_mpc_best.jls.

## Confirmation 1 — forecast-value curve STILL FLAT (test/forecast_value.jl)
    fcst_skill 0.00 0.25 0.50 0.75 1.00  ->  2.83 2.87 2.85 2.84 2.79
Even the distilled agent does NOT respond to forecast accuracy. The anticipation it
learned is NOT weather-forecast-driven -- it is a better SEASONAL + STATE-conditioned
policy (day-of-season, crop stage, current climate); the forecast-specific deviations
in the expert's actions are weak/noisy and averaged out in BC. So it learned the
CLIMATOLOGICALLY OPTIMAL PROGRAM SHAPE, which alone gives 78% of the ceiling.
=> short-term forecast has little marginal value in light-limited autumn; DO NOT
extend _FC_H 3->7.

## Confirmation 2 — vs 6 real AGC-2018 growers (test/compare_vs_agc.jl)
Real growers' data, OUR economics, density 2.5, elec blended 0.25:
    grower/agent   kg/m2  heat  elec  CO2   netEUR/m2(@0.25)
    AiCU           37.3  116.8 118.3  9.9   -15.19
    Croperators    48.0  233.9 183.3 14.2   -33.75
    DeepGreens     35.8  486.9 157.0 16.7   -61.59
    Reference      49.4  158.2 149.3  9.8   -15.88
    Sonoma         51.3  127.9 184.2 10.3   -20.41
    iGrow          46.8  137.3 145.9  9.2   -15.33
    AGENT ref-crop 28.3   84.7  29.5  5.0    +4.03
    AGENT native   22.8   63.2  31.6  1.2    +1.64
Agent grows LESS but spends far less. Not proof it grows better.

## Confirmation 3 — electricity-price sweep (the honest verdict)
net profit vs elec price [EUR/kWh], selected rows:
    price          0.05   0.10   0.15   0.20   0.25   0.30
    Reference     13.98   6.52  -0.95  -8.42 -15.88 -23.35
    Sonoma        16.43   7.22  -1.99 -11.20 -20.41 -29.62
    iGrow         13.85   6.55  -0.74  -8.04 -15.33 -22.63
    AGENT ref      9.94   8.46   6.99   5.51   4.03   2.56
    AGENT native   7.96   6.38   4.80   3.22   1.64   0.06
Agent profit is nearly FLAT vs price (barely uses power); growers steeply
price-sensitive. CROSSOVER ~0.10-0.12 EUR/kWh. At 2018-era cheap power (~0.05) the
BEST growers (Sonoma 16.4, Reference 14.0, iGrow 13.9) OUT-EARN the agent (9.9).
Conclusion: the agent is PRICE-ROBUST, not a better grower; it wins only where power
is expensive, by declining lighting that doesn't pay in a light-limited autumn.
Caveat: agent policy is FIXED (trained at ~0.25); a fair cheap-power test needs a
retrain at that price (it would then light up).

## Confirmation 4 — OOD diagnostic (test/plot_actions.jl [actor] [year] [native|matched])
Matched-crop trajectory: agent stable ~80 d, then erratic in the final weeks --
night setpoint jumps ABOVE day setpoint, heating spike, air temp crash to ~8 C,
yield flatlines. NATIVE-crop run under the SAME weather year is smooth to the end
(temp held ~17 C through the same day-82 outside cold dip to ~5 C). => the anomaly
is OUT-OF-DISTRIBUTION (crop-state obs, esp. cumulative yield, exceed the agent's
training range on the higher-yield Reference crop and drive policy extrapolation),
NOT a twin/weather bug. Agent reliable in-distribution; matched-crop 28 kg/m2
UNDERSTATES its control quality. Fix = short fine-tune on the matched crop.

## RL line — conclusion (paused here)
CEM 0.12 -> reactive 1.66 -> MPC-distilled 2.83 (climatological optimum) -> ceiling
3.65. Forecast responsiveness ~nil. Price-robust, not better-grower. Reliable
in-distribution. Consolidated report written (docs/greenhouse_control_report.pdf,
delivered) with all four confirmations + price sweep + OOD diagnostic; To-Do led by
PRICE-ADAPTIVE agents (prices in the observation, trained across price regimes).

## Next
git monorepo reorg (GreenhouseCore + Inference + Simulator + Control) once user
provides GitHub credentials. Then goal #1: t-walk inference + probabilistic 3-5 wk
production forecast.
