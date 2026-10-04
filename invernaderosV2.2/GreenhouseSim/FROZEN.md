# GreenhouseSim V2.2 — FROZEN reference model (project goal #1)

This package is the validated **model engine**: greenhouse climate physics +
FvCB photosynthesis + Marcelis-style cucumber growth, plus the calibration
machinery. It is **frozen** — the control/agent work (goal #2) lives in a
separate package (`invernaderosV3.0/GreenhouseControl`) that reuses this engine
by `include`, without copying it. Do not add agent/control code here.

## What it provides

- Fast, type-stable forward map: `climate_rhs`, `assimilation`, `rhs!`,
  `forward_map` (physical constants in one concrete NamedTuple; no globals).
- Climate-parameter inference (t-walk): `test/run_inference.jl`.
- Cucumber growth model: `grow!`, `simulate_growth`, `simulate_growth_measured`.
- Growth calibration on measured AGC-2018 data: `test/calibrate_growth.jl`,
  `test/calibrate_all_teams.jl`; fast visual `test/check_growth.jl`.

## Calibrated / fixed parameter status (cucumber, AiCU reference)

Identified from data:
- `J_max ≈ 1.15e-4` (≈115 µmol; light-limited regime), `node_rate ≈ 0.091`,
  `veg_sink_max ≈ 21`, `set_start_day ≈ 16`, `LAI_max ≈ 2.5`.

Fixed at physical / literature values (see configfiles/cucumber_growth.json):
- `k_ext = 0.75` (canopy extinction), `asrq = 1.4`, `SLA = 0.03`,
  `V_cmax25` (non-identifiable under heavy CO2 dosing — crop is light-limited),
  `top_day = 104`, `set_r_low/high` (flush shape), `stem_density = 2.5`.

Weakly identified (expected): `set_rate`, `Wf_max` (yield can't separate
many-small vs few-large fruit).

Cross-team validation: all 6 AGC growers fit with the same model/priors, yields
within ~10%; development params cluster; source params vary ~2:1 (partly the
plant-density confound). See project notes 05 and 06.

## Backward-compatible additions since first freeze

- `stem_density` growth parameter (default 2.5 = REF_DENSITY ⇒ scaling factor 1,
  so all prior results reproduce exactly). Scales the per-area capacities
  (`veg_sink_max`, `set_rate`, `LAI_max`) by `stem_density / 2.5`. Enables
  variable-density studies and de-confounds cross-team calibration when real
  per-team densities are supplied (to BOTH `gp.stem_density` and
  `load_cropmanagement`'s `stem_density`).

## Physics corrections (controlled, calibration-safe)

- **Ventilation rate sign (climate.jl `f7`)**: removed a `* sign(T2 - I5)` factor
  that made the roof-ventilation RATE negative when inside was colder than outside.
  A ventilation rate is a non-negative air-exchange rate; transport direction is
  carried by the (T2-I5)/(C1-I10)/vapour differences in the balances. The old sign
  turned the CO2 and vapour ventilation terms into a positive feedback (blow-up)
  at realistic leakage. **Identical whenever T2 >= I5** (all venting in the summer
  calibration), so calibrated forward map / parameters are unchanged — verify with
  `test/benchmark.jl` (regression) after this change.

- **Roof-vent sensible heat (climate.jl `h7`)**: the air-temperature ventilation
  loss `h7` was `(f2 + f3 + 0.5·f6)` — side vent + forced + half-leakage — and
  omitted `f4` (the roof/buoyancy vent), even though the CO2 (`o5`) and vapour
  (`p5`) balances both include `f4`. A ventilation air-exchange rate carries
  sensible heat whatever drives it, so the roof term belongs here too. Fixed to
  `(f2 + f4 + f3 + 0.5·f6)`. This is the root cause of the closed-loop twin's
  ~40 °C excursions on **calm** days (side vent `f5 ∝ windspeed → 0`, so with no
  roof-vent cooling the greenhouse could not shed the solar load). **Unlike the
  sign fix, this DOES change the forward map whenever the roof vent is open
  (`U8 > 0`, i.e. most summer days), so it is NOT calibration-neutral**: re-run
  `test/benchmark.jl` to quantify the shift and re-run climate inference
  (`test/run_inference.jl`). Expected consequence: the leakage `nu4/nu4CO2`,
  previously inflated ~100× above Vanthoor to force daytime cooling the model
  otherwise lacked, should relax toward physical values (~1e-4).

  **VALIDATED (done).** The twin leakage sweep (V3.0 `sweep_leakage.jl`) shows the
  coupled twin is stable and better-behaved (CO2 holds setpoint, no night undershoot)
  all the way down to 1e-4. A t-walk MCMC on a real wide-open venting day under the
  corrected physics (`test/run_inference_venting.jl`, AGC-2018 Reference/Growers,
  Tair+RH, daytime hours) gives a `nu4` posterior concentrated at ~2-8e-4 with the
  old 1e-2 fully excluded. The JSON default `nu4 = nu4CO2` is therefore restored to
  **1e-4** (was 1e-2; backup `constants_climate.json.bak_preleak`). NOTE: this changes
  the default-parameter forward map on leakage-sensitive windows, so re-baseline
  `test/benchmark.jl` (its old reference assumed 1e-2). `gamma4` leaned to its lower
  bound in the same fit (weakly identified / partly absorbing the daytime thermal-mass
  amplitude gap) and was left at its literature value.

## Floor/soil thermal-mass state (added)

The 4-state climate model treated soil as a fixed 18 C input, so it had no heat
store: the air over-swung (too hot at midday, too cold at night) and the fit pushed
`gamma4` to its bound to fake the missing damping. Added ONE dynamic floor state
`T5` (appended LAST to each RHS u-vector; `climate_rhs` itself unchanged): the floor
is fed to the existing soil couplings in place of the fixed input, and a single
conductance `k_ground` charges it from the air and discharges it to the deep-soil
anchor (fixed capacity `C_soil = 1.5e5 J/m2/K`). The air loses exactly what the floor
gains (`floor_balance` in forward_map.jl; used by `rhs!`, `rhs_dyn!`, `rhs_control!`).

Calibrated on the AGC venting day (`test/run_inference_venting.jl`, daytime hours):
**`k_ground = 4.33 W/m2/K`, 95% CI [1.45, 7.82]** -- cleanly identified (an earlier
`h_soil`+`C_soil` pair was degenerate and railed; collapsed to this one conductance).
The daytime air-temperature peak drops from ~38 C to the measured ~28 C. `nu4` stays
`1.65e-4` (leakage conclusion intact); `gamma4` still ~24 (a genuine vapour-balance
preference, not thermal-mass compensation -- it did not move when the floor was added).
LIMITATION: daytime-only data pin the air<->ground conductance but not the storage
capacity, so `C_soil` is FIXED at a physical value rather than inferred.

## Cloud-derived sky temperature (added, opt-in)

`load_weather(...; sky_from_clouds=true)` computes the effective sky (radiative)
temperature per hour from air temperature, RH and total cloud cover (meteo col 9)
via a Berdahl-Martin clear-sky emissivity raised toward 1 (Tsky -> Tout) as cloud
increases -- instead of the fixed -0.4 C clear sky. Default is `false` (exact old
behavior; the frozen inference paths are unchanged).

WHY: a season heat-balance audit of the control twin (V3.0 `test/heat_audit.jl`)
showed FIR-to-sky was the dominant loss (394 kWh/m2, ~164 W/m2 avg) and the pipe
heating (208 kWh/m2) was mostly replacing heat radiated to a permanently clear cold
sky. With the cloud-derived sky, FIR drops to 230 and heating to 129 kWh/m2
(1.29 kWh/m2/day vs the AGC-2018 Reference grower's 1.18) -- the ~1.8x over-heating
is resolved. The twin/RL runners (V3.0) pass `sky_from_clouds=true`. Consistent with
the venting calibration, which already used a measured (pyrgeometer) sky.

## Energy-screen emissivity (added)

The thermal-screen sky-FIR term `r12 = epsil5*(U1*tau2)*SkyFIR` used epsil5=1, tau2=1
-- a FULL-emissivity shade screen that *radiates more* to the sky when deployed, so it
retained no heat (a season audit found r12 = 152 kWh/m2, the dominant FIR loss, growing
when the screen closed). Real energy screens are aluminized/low-e. Added a separate
`eps_screen = 0.4` used only in r12 (screen deployed, U1>0), so closing the screen now
cuts radiative loss. Calibration-safe: U1=0 on the venting-calibration day, so r12=0 and
nothing calibrated moves. Effect: with the V3.0 night energy screen, r12 152->61,
heating 140->100 kWh/m2. `eps_screen` is calibratable.

## Known limitations (documented, not blocking)

- Absolute source/allocation params ride on the assumed 2.5 stems/m² until real
  densities are used.
- Weekly harvest-flush amplitude is approximate; `V_cmax25` needs Rubisco-limiting
  data to identify.
