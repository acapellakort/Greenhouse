# Goal #2 — control agent: architecture & versioning

## Versioning decision

- **V2.2 (`invernaderosV2.2/GreenhouseSim`) is FROZEN** = the validated model engine
  (physics + FvCB + cucumber growth + calibration). Goal #1 deliverable. See its
  `FROZEN.md`. Do not add control/agent code here.
- **V3.0 (`invernaderosV3.0/GreenhouseControl`) is NEW** = the control/agent layer
  (goal #2). It **reuses V2.2 by `include`, never copies** — single source of truth,
  avoiding the 2.0/V2.1/3.0 drift. Reuses 2.0's proven Julia env.
- **claude.ai project unchanged** — "Greenhouses" spans both goals.

## V2.2 backward-compatible change before freeze: variable density

`stem_density` growth param (default 2.5 = REF_DENSITY ⇒ factor 1). Scales per-area
capacities (`veg_sink_max`, `set_rate`, `LAI_max`). Enables density studies + de-confounds
cross-team calibration when real densities are supplied.

## Control architecture (hierarchical) — validated against AGC-2024 doc (note 08)

**Agent acts in SETPOINT space; a fixed PID layer executes.** PI heating +
proportional-band venting (dead-zone offset) + PI CO2 + rule screen/lights, on **air temp T2**. The
AGC-2024 simulator uses the identical scheme — independent confirmation.

## Stage 1 (DONE + VALIDATED) — control layer + closed-loop RHS

`src/climate_computer.jl`: `Setpoints` + `daynight_setpoints`, `vent_control`, `PIGains`,
`rhs_control!` (state `[T1,T2,C1,V1,tempInteg,co2Integ]`; continuous leaky-integral PI for
heating → pipe temp I3 and CO2 → valve U10; vents; screen; lights), `run_control`.
Demo `test/run_control.jl` → tracking + actuator plot.

**Result:** the twin tracks a programmed climate cleanly — air holds day/night Tset
(~18/21 °C), CO2 enriches to setpoint, actuators modulate. Gains `gT=(0.6,0.02,5e-4)`,
`gC=(0.01,2e-4,5e-4)` work well.

### KEY FIX found here (root cause of both tracking failures)

The heat/CO2 **leakage coefficients `nu4`, `nu4CO2` were inflated ~100×** (0.01 vs
Vanthoor ~1e-4) in the 2.0 model *to stabilize the ODE* (per its README). That artificial
leak caused BOTH the night-temperature undershoot AND CO2 that wouldn't enrich. Overriding
them **10× lower (1e-3)** in the V3.0 twin (`update_params`, V2.2 stays frozen) fixed both
at once. Also `psi2` bumped to 27800 (AGC CO2 capacity 100 kg/ha/h). If pushed toward full
physical 1e-4 and the ODE stiffens, switch to a stiff solver (TRBDF2/Rodas).

## Remaining minor items (control side)

- Daytime ~1–2 °C overshoot on hot midday (physical vent limit). AGC refinements would
  smooth it: P-band widens with cold/wind, day/night vent offset, radiationInfluence
  (raise Tset with light), setpoint slope limits (2 °C/h, 500 ppm/h). See note 08.
- Energy screen (AGC has a dedicated 70%-transmission night screen) — our single screen
  U1 is cruder; optional refinement.

## Next

- **Stage 2**: couple growth into the loop (daily: PID climate day → daily assimilate →
  `grow!` → LAI feeds back) = full closed-loop twin (setpoints → climate → growth → yield).
- **Stage 3**: RL environment — action = setpoints, reward = net profit from the AGC-2024
  cost coefficients (note 08), with cucumber economics. Then the agent.

## Run

    cd invernaderosV3.0/GreenhouseControl
    julia --project=. test/run_control.jl
