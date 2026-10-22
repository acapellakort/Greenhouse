# GreenhouseTwin — the n79 model without flags

Climate (Vanthoor-type ODEs) + FvCB photosynthesis + Aaron Vélez's generic cucumber/canopy crop model for the AiCU
compartment of the Autonomous Greenhouse Challenge 2018. It is the V2.2 engine (tag `v2.2-flags-final`, commit 614a178)
with the **n79 configuration resolved into one equation set**: every process listed below is always on; there are no
behaviour flags. Equations: `docs/greenhouse_model_n79.tex`.

## What "n79" resolves to (V2.2 flag → what the code now always does)

| V2.2 setting | in GreenhouseTwin |
|---|---|
| use_cover_node, cover_area 1.6, cover_cond | algebraic cover node with condensation and its latent heat; air–cover term × 1.6 |
| solar_cons, solar_mass_share, par_below (SOLARC=1) | energy-conserving solar column (roof → energy screen → blackout screen → canopy → floor); mass takes f_mass of floor-level solar; plant light = PAR on both canopy faces + lamp PAR |
| phys_light, photo_lai, rd_q10 | photosynthesis on plant PAR × effective leaf area (1−e^(−0.7L))/0.7, leaf dark respiration Q10 2 on canopy T |
| stomata_fix, dew_bl | stomatal resistance with day/night coefficients in the right order and CO₂ in ppm; dew sees only the boundary layer |
| vapour_fix, cond_heat (f_cond 0) | screen condensation × U1 with the right temperature difference; its latent heat not returned to the air |
| vent_split, vent_screen (open loop) | logged leeward + windward roof vents (one equivalent opening), throttled by a closed energy screen |
| leak_stack, leak_shared | leakage ν4·√(max(v,0.25)² + c_s·|ΔT|), full leakage in heat, vapour and CO₂ |
| heat_ctrl, hc_log_pipe | logged pipe where it is warm, PI controller to the heating setpoint elsewhere |
| floor_full, solar_split, deep_k, soil_dyn | closed floor balance (convection, conduction, canopy and cover long-wave, solar), own floor→deep conductance, deep-soil state |
| mass_node, crop_cap | internal mass state; its capacity includes the standing crop fresh mass |
| lamp_full (117 W m⁻²) | full lamp power: canopy share lamp_fcan·(1−e^(−β1 L)), rest to the air |
| co2_resp | crop (non-leaf) respiration as a CO₂ source, Q10 2 |
| air_cap_mult | effective air heat-capacity multiplier (inferred) |
| fixed geometry (Note 78) | φ2 5.4 m, vent area 0.30, pipe 1.0 m m⁻² (one circuit), PAR transmission 0.65, CO₂ valve 12.8 g m⁻² h⁻¹ |
| removed | vent controller, pipe lag, blackout insulation, s_rad, dday, hc_dT, lamp screen loss, V2.0 lamp terms, V2.0 solar terms, legacy crops |
| kept as inputs (zero in AiCU) | fan-pad, mechanical cooling, direct heater, fog, side vents, shading, forced ventilation |

State: `u = [T1 canopy, T2 air, C1 CO₂, V1 vapour pressure, e_int, T5 floor, Td deep soil, Tm internal mass]`.
Crop: offline coupling — calibrated on measured climate; `test/crop_series.jl` writes the LAI / respiration / fresh-mass
series the climate model reads.

## Files

| path | content |
|---|---|
| `src/climate.jl` | `model_terms` (every flux and derivative), `rhs!`, cover node |
| `src/solar.jl` | solar column |
| `src/photosynthesis.jl` | FvCB canopy assimilation |
| `src/crop.jl`, `src/crop_measured.jl` | Aaron's crop + LAI model; measured-climate driver, observations, likelihood, climate series |
| `src/agc_data.jl` | AGC-2018 inputs and observations (n79 conventions) |
| `src/climate_calibration.jl` | blocks, likelihood, season evaluation |
| `configfiles/climate_n79.json`, `crop_n105.json` | every parameter: value, units, group, fixed/inferred, prior, source |
| `tools/make_params.py` | how the JSONs were built from V2.2 (provenance) |
| `tools/audit.py` | static checks: no flags left; JSON ⇔ code parameters |

## Commands (from this folder)

First time: `julia --project=. -e 'using Pkg; Pkg.instantiate()'` (Manifest copied from V2.2).

```
julia --project=. test/equivalence_synthetic.jl            # V2.2-n79 vs this code, synthetic inputs (~1 min)
julia --project=. test/equivalence_n79.jl                  # V2.2-n79 vs this code, AGC data (minutes)
julia --project=. test/check_solar_column.jl
julia --project=. test/benchmark.jl
MODE=centre julia --project=. test/calibrate_crop.jl
SEEDS=2 BURN=150000 julia --project=. test/crop_series.jl  # crop series from the V2.2 n105 chains
CROPTAG=n105 RELTOL=1e-4 julia --project=. test/calibrate_climate.jl smoke
python3 tools/audit.py
```

## Verification status (21 Oct 2026; cloud sandbox, Julia 1.10.4 — to be repeated on the Mac)

| check | result |
|---|---|
| `tools/audit.py` | PASS |
| `test/equivalence_synthetic.jl` | 35 / 35 |
| `tools/check_tex.py` (document vs code, 60 terms) | 60 / 60 |
| `test/equivalence_n79.jl` (AGC data) | see Note 109 |
| `test/check_solar_column.jl` | 39 / 39 |

**Numerical noise.** Trajectories and the likelihood react to rounding-level differences (discontinuous inputs and
switches). At reltol 1e-4 the energy moves by up to ~60 for a 1e-14 parameter change; at 1e-6 by ~1; at 1e-7 by ~0.1.
This is inherited from V2.2-n79. Deal with it before long climate MCMC runs (Note 109).
