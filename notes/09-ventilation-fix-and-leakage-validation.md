# 09 — Ventilation physics correction and leakage re-validation

This note records a two-part correction to the V2.2 climate physics (both inherited
from the original 2.0 model) and the data-based re-validation of the heat-leakage
coefficient that followed. Context: building the Stage-2 closed-loop twin
(V3.0 GreenhouseControl) surfaced ~40 °C air-temperature excursions that were not
physical.

## The two ventilation bugs (climate.jl)

1. **Ventilation-rate sign.** The roof-ventilation rate `f7` carried a `* sign(T2 - I5)`
   factor that made the air-exchange RATE negative when inside was colder than outside.
   A ventilation rate is a non-negative buoyancy magnitude; transport direction is
   already carried by the (T2−I5)/(C1−I10)/vapour differences. At realistic leakage the
   sign turned the CO2 and vapour vent terms into positive feedback → CO2 blow-up.
   Removed. Identical whenever T2 ≥ I5 (all summer venting), so calibration-neutral.

2. **Roof vent missing from the air heat balance.** The sensible-heat vent loss
   `h7 = ρ·α·(f2 + f3 + 0.5·f6)·(T2−I5)` included the *side* vent `f2` (wind-driven,
   `f5 ∝ windspeed`) but omitted `f4`, the *roof/buoyancy* vent — even though the CO2
   (`o5`) and vapour (`p5`) balances both include `f4`. So the roof vent flushed CO2 and
   moisture but removed no heat. On calm days the side vent → 0 and the greenhouse could
   not shed its solar load → runaway to ~40 °C. The spikes correlated with **low wind,
   not high outside temperature** (day-8 case: 41 °C inside vs 22.7 °C outside max, on the
   calmest day). Fixed to `h7 = ρ·α·(f2 + f4 + f3 + 0.5·f6)·(T2−I5)`. Confirmed the same
   omission exists in the original 2.0 `climate_functionsV3_NoUnits.jl` (line 145), so the
   V2.2 port was faithful — this is an original-model bug. Adds a stabilizing cooling term
   (cannot reintroduce blow-up). Unlike the sign fix, it DOES change the venting-day map.

## Leakage restored to physical (nu4)

The JSON comment on `nu4` stated it plainly: Vanthoor's value is 1e-4, "lo hacemos mas
grande para estabilizar el sistema de ODE" — inflated to 1e-2 only to stabilize the ODE,
i.e. to paper over the bugs above.

- **Twin leakage sweep** (`V3.0 test/sweep_leakage.jl`): with both fixes in, the coupled
  twin is stable from 1e-3 down to the physical 1e-4. As leakage drops, mean CO2 rises
  toward the 1200 setpoint and the night temperature stops undershooting — both improve.
- **Forward validation** (`V3.0 test/validate_leakage_agc.jl`): driving the corrected
  model with AGC-2018 measured weather + measured roof apertures, the daytime venting
  temperature tracks the measured Tair and the 1e-4 / 1e-3 / 1e-2 curves are
  indistinguishable → leakage buys no fidelity.
- **Venting-day MCMC** (`V2.2 test/run_inference_venting.jl`): t-walk on a real wide-open
  venting day (AGC Reference/Growers, Aug 21), scoring Tair + RH over daytime hours only.
  `nu4` posterior at ~1.6e-4 with the old 1e-2 fully excluded.

**Decision:** JSON default `nu4 = nu4CO2` restored to **1e-4** (backup
`constants_climate.json.bak_preleak`). Re-baseline `test/benchmark.jl` — its old
reference assumed 1e-2. Growth calibration is unaffected (uses measured climate, not the
ODE). New AGC climate loaders (weather from meteo.csv timestamped by GHtime; controls
mapping VentLee/Ventwind → roof U8 with side U6 = 0; Tsky from measured net longwave
Pyrgeo) live in the test scripts and are reusable.

## Floor/soil thermal-mass state (added, calibrated)

The 4-state model treated soil as a fixed 18 °C input, so it had no heat store: the air
over-swung (too hot at midday, too cold at night) and the MCMC pushed `gamma4` to its
bound to fake the damping. Added ONE dynamic floor state `T5`, appended LAST to each RHS
u-vector so no existing index moves; **`climate_rhs` is byte-identical** — the floor is
simply fed to the existing soil couplings in place of the fixed input. A shared helper
`floor_balance` (forward_map.jl, used by `rhs!`, `rhs_dyn!`, `rhs_control!`) adds one
conservative convective exchange (air loses exactly what the floor gains) plus a deep-soil
anchor.

First attempt inferred `h_soil` + `C_soil`; they were **degenerate** (h_soil railed high,
C_soil railed low — daytime-only data pin the air↔ground coupling but not the storage
capacity). Reparametrized to a **single conductance** `k_ground` with a fixed physical
capacity `C_soil = 1.5e5 J/m²/K`. Result (venting-day MCMC, daytime hours):

- **`k_ground = 4.33 W/m²/K`, 95% CI [1.45, 7.82]** — cleanly identified, no railing.
- Daytime air-temperature peak drops from ~38 °C to the measured ~28 °C; posterior-median
  model sits on the measured line through the scored window.
- `nu4` stays **1.65e-4** (leakage conclusion intact); `gamma4` still ~24 — it did NOT move
  when the floor was added, so that low value is a genuine vapour-balance preference
  (well-watered cucumber transpires heavily), not thermal-mass compensation.

JSON default `k_ground = 4.33`. LIMITATION: daytime-only data cannot identify the storage
capacity, so `C_soil` is fixed at a physical value rather than inferred; pinning it would
need a night-inclusive fit, which needs better night heating data than the sparse pipe log.

## Status

Goal #1 is complete: the inference twin is fast (V2.2 type-stable), growth-calibrated,
physics-corrected (both vent bugs), leakage-validated at the physical value, and now
thermal-mass buffered with a single calibrated conductance. Ready for Stage 3 (RL
net-profit environment).
