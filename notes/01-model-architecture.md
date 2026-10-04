# Greenhouse Model — Architecture & How to Run

_Notes folder for the Greenhouses project. This is the first note: a map of the
existing Julia software (climate physics + photosynthesis), how each version runs,
and where the plant/vegetable growth model must plug in._

Source folder on santiagos-laptop-local:
`~/Documents/trabajo/academia/unam/proyectos investigacion/bio/invernaderos`

---

## 1. The shared model (same physics in every version)

A **Vanthoor-type greenhouse climate model coupled to an FvCB
(Farquhar–von Caemmerer–Berry) C3 photosynthesis model**, solved as a system of ODEs.

**State vector** `u = [T1, T2, C1, V1, …]`:

| Symbol | Meaning | Units |
|--------|---------|-------|
| T1 | Canopy temperature | K |
| T2 | Greenhouse air temperature | K |
| C1 | CO₂ concentration | mg/m³ |
| V1 | Air vapour pressure | Pa |
| … | PID integral-error terms + stored assimilation `A` | — |

**Climate RHS** (`rhs_fast` in 2.0/V2.1, `climate_model_rhs` in 3.0) — four coupled balances:

- `dT1` — canopy energy balance (PAR/NIR/FIR radiation, sensible heat to air, latent heat of transpiration, pipe radiation, sky/floor FIR).
- `dT2` — air energy balance (canopy/pipe/pad/mech-cool/floor exchange, ventilation losses, lamps, global radiation on structure).
- `dC1` — CO₂ mass balance (dosing, pad-fan, canopy uptake `A`, ventilation/leakage).
- `dV1` — water-vapour balance (transpiration, pad, fog, blower, ventilation, condensation terms).

All physical constants (α, β, γ, ν, η, τ, φ, ψ families) load from JSON:
`climate_parameters.json` / `constants_climate.json`.

**Photosynthesis** (`compute_assimilates` / `AcropFast`) — full C3 FvCB model.
Returns assimilation `A` (mg CO₂ / m² / s) as `LAI × min(Ar_j, Ar_c, Ar_p)`:

- `Ar_c` — RuBisCO-limited
- `Ar_j` — electron-transport (light) limited
- `Ar_p` — product (triose-phosphate) limited

Temperature-dependent `V_cmax`, `J_max`, `K_C`, `K_O`, CO₂ compensation point `Γ*`.
Params from `photosynthesis_parameters.json` / `constants_cropphoto.json`.
Stomatal conductance is currently a **constant** `g_eff()` (Bush 2023); an
environment-dependent `r_s(T, I, C, V, LAI)` exists but is **not wired in**
(README: "FALTA: modelo con dependencia en variables ambientales").

**Coupling & structure:**
- `A` feeds back only into `dC1` (canopy CO₂ uptake). One-directional.
- **LAI is an input, not a state** — hardcoded `0.5` in 3.0, or from weather interp `I1` in 2.0.
- **Daily outer loop**: each day solved with `BS3()`/`Tsit5()`, end-state carried to next day.
  The README notes this was designed so a growth model can act **once per day** in that loop.
- Controls `U1…U12`: U1 thermal screen, U2 pad-fan, U3 mech cooling, U4 air heater,
  U5 shade screen, U6 side vents, U7 forced vent, U8 roof vents, U9 fog, U10 CO₂ source,
  U11 pipe setpoint, U12 lights. Driven by PID + `climateComputerOrders.json` in simulation,
  or by measured control CSVs (interpolated) in inference.
- Solver note (README): `Tsit5()` chosen for robustness; `Trapezoid/QNDF/FBDF` fail to converge.

---

## 2. The three versions

### invernaderos2.0 — mature, working version (the workhorse)
- `main_program.jl` — forward simulation with live PID / rule-based controls.
- `for_inference_V4.jl` — **inverse problem**: reads measured controls from CSV (as
  time-driven inputs), simulates forward, runs **t-walk MCMC (`JTwalk`)** to infer
  parameters (default `beta2, gamma3, gamma4, nu4`) against observed T1/T2/RH/CO₂;
  saves chain + posterior figures to `figures/`.
- Parameters are **global constants** set in `utils/import_constants_NoUnits.jl`.
- This is the version to trust for results today.

### invernaderosV2.1 — clean refactor of 2.0, WORK IN PROGRESS (won't run as-is)
- Proper Julia package `GreenhouseSim` (`Project.toml`, module, `ModelParams` struct,
  `Controls` module, `SimulationContext` — no globals). Good target architecture.
- **Blockers before it runs:**
  - `run_inference.jl` calls `GreenhouseSim.load_weather_data(...)` — **does not exist** in the package.
  - Expected `utils/weather_data.csv` and `data/observed_data.csv` are **not present**.
  - `FM_SimLoop!` has **no return statement** (code says "Add return logic"), so `energy()` fails.
  - `test/` is empty.
- Intended entry: `run_inference.jl`.

### GreenHouse_3.0 — GUI for direct simulation
- Cleanest module layout: `climateEquationsRHS.jl`, `photosynthesisEquations.jl`,
  `controlFunctions.jl`, `parameterLoader.jl`.
- `GUI.jl` — **Dash web app** (`HarverstSim` module): choose start date + number of days,
  run simulation, animated overlay on a greenhouse image + time-series plots.
- `simulation_main.jl` — headless equivalent that produces plots.
- ⚠️ The `src/` subfolder is an **older stale copy** of the same code — ignore for edits.

---

## 3. How to run each

| Task | Commands |
|------|----------|
| 3.0 GUI | `cd GreenHouse_3.0/code` → `julia --project=. -e 'using Pkg; Pkg.instantiate()'` → `julia --project=. GUI.jl` → open `http://localhost:8050` |
| 3.0 headless plots | `julia --project=. simulation_main.jl` |
| 2.0 forward sim | `cd invernaderos2.0/code` → `julia --project=. main_program.jl` |
| 2.0 inference (t-walk) | `julia --project=. for_inference_V4.jl` (needs `../observed_data.csv`) |
| V2.1 | intended `cd invernaderosV2.1/GreenhouseSim` → `julia --project=. run_inference.jl` — **blocked** until weather loader + data files + `FM_SimLoop!` return are added |

---

## 4. The gap — the growth model (project goal)

There is **no plant/vegetable growth model yet**. Photosynthesis produces `A`, but:
- `A` only feeds the CO₂ balance — no biomass accumulates.
- LAI is a fixed input, not a growing state.

**Insertion point** — the daily outer loop:

```
integrate A over the day  →  daily assimilate gain
      → partition dry matter (leaves / stems / fruit)
      → update biomass + LAI
      → feed new LAI into next day's climate + photosynthesis
```

That closes the loop into a digital twin.

---

## 5. Open decisions

- **Where to build the growth model:** onto 2.0 (proven, but global-heavy) vs. first
  finish the V2.1 package (cleaner foundation) and build there.
- Whether to wire in the environment-dependent stomatal resistance `r_s` (already coded, unused).
- Growth-model formulation to adopt (e.g. TOMGRO / reduced-state dry-matter partitioning).

_Next notes to add here: full parameter inventory, RHS equation write-up, V2.1 fix list, growth-model design._
