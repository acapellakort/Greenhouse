# Growth-parameter calibration (t-walk on measured cucumber yield)

Fits the cucumber growth parameters to the real AGC-2018 yield, using the same
t-walk MCMC as the climate inference. See notes/04 for the growth model.

## Design: drive photosynthesis from MEASURED climate; multi-objective likelihood

Daily gross assimilate is computed directly from the measured greenhouse climate
(`Tair`, `CO2air`, light) via FvCB `assimilation` (no ODE; ~30k evals/season,
<1 ms). We fit THREE observables jointly, each with its own inferred sigma:
1. weekly incremental yield (Production.csv `Total_Prod_cum`);
2. cumulative leaf/node number (CropManagement.csv `N_leaves`) → development;
3. derived standing LAI (from N_leaves) → vegetative allocation.

Weekly increments (not cumulative) avoid autocorrelation-inflated posteriors.
LAI derivation assumes leaf_area 0.05 m², stem_density 2.5/m² (ASSUMED), N_window
20 → plateau ~2.5 (`load_cropmanagement`; week→day vs day0 = 2018-08-14).

## Iteration history (AiCU)

- **Iter 1–2** cumulative/weekly, QoI incl. `asrq` → overshoot ~50%; `asrq` railed.
  Cause: canopy photosynthesis scaled LINEARLY with LAI.
- **Iter 3** Beer's-law `lai_eff=(1−e^(−kL))/k` + inferred `set_start_day`.
  Overshoot →~16%; onset aligned; `asrq` still railed at 1.9.
- **Iter 4** fix `asrq=1.4`, infer `k_ext`. Cumulative fit good but `k_ext~1.0`
  and `set_rate/veg_sink/Wf_max` uniform (yield alone can't identify partitioning).
- **Iter 5** MULTI-OBJECTIVE (add leaf + LAI). `veg_sink_max` peaked ~21,
  `node_rate` sharp ~0.091, leaf/LAI good. But `k_ext` railed at 1.3 ⇒ leaf-level
  FvCB source ~30–40% too high.
- **Iter 6** fix `k_ext=0.75`; infer `V_cmax25`, `J_max`. `J_max` identified
  **~1.15e-4 (≈115 µmol)**; `V_cmax25` **flat/unidentified** — crop is light-limited
  throughout (heavy CO2 dosing lifts the Rubisco ceiling), so only J_max binds.
- **Iter 7** infer only `J_max`; fix `V_cmax25` at literature (non-binding). Good
  fit at physical params (cum ~35 vs 37).
- **Iter 8 (current) — STRUCTURAL TIGHTENING** (see notes/04 growth.jl):
  - **Topping** (`top_day=104`): node/leaf formation + fruit set stop late season.
    Fixes leaf plateau (model now plateaus ~106 at day 104, matches obs ~105) and
    bends the late-season yield tail.
  - **Carbohydrate-regulated fruit set** (`set_r_low=0.6`, `set_r_high=1.4`):
    new set is suppressed when supply/demand is low and resumes after harvest;
    with the ~2-week fruit-growth delay this produces **harvest flushes** (model
    weekly yield now oscillates ~1.5–2 wk, amplitude ~0.8, realistic).
  - Flush params are FIXED (not inferred) — the weekly point-wise SSE is phase-blind
    and would just suppress the oscillation. `check_growth.jl` = fast 1-season
    visual (no MCMC) to tune them.
  - Side effect: set-suppression removes some fruit ⇒ level dropped (~31 at OLD
    point estimates). Re-run calibration so `set_rate` (and `J_max`) re-settle and
    restore ~37.

## Parameter status

- Identified: `J_max` (~115 µmol), `node_rate` (~0.091), `veg_sink_max` (~21),
  `set_start_day` (~16), `LAI_max` (~2.5, rides on assumed density).
- Fixed: `k_ext=0.75`, `asrq=1.4`, `SLA=0.03`, `V_cmax25` (non-ID),
  `top_day=104`, `set_r_low/high` (flush shape).
- Weakly identified (expected): `set_rate`, `Wf_max`.

## Scientific finding (relevant to control agent, goal #2)

Under AiCU's CO2 strategy the crop is **light-limited**, so yield responds to
LIGHT (lamp use, canopy light capture, `J_max`), NOT to Rubisco capacity.

## Scripts

- `test/check_growth.jl` — fast 1-season visual at point estimates (tune structure).
- `test/calibrate_growth.jl` — full multi-objective t-walk.

## Next

- Re-run calibration on the tightened model (level back to ~37).
- Calibrate other teams (parameter spread across growers) — confirm reasonable.
- Then: fully-coupled ODE twin (AGC meteo + mapped controls) for the RL phase (goal #2).

## Ops note

- `device_commit_files` sometimes reports success but does NOT write; verify
  on-disk with grep after each push and retry if stale.
