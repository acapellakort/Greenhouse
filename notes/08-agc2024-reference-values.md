# AGC-2024 simulator reference values (equipment, control, economics)

Source: "Part B — Climate and crop control Challenge (Simulator)", 4th Autonomous
Greenhouse Challenge, WUR, 2024-03-15 (dwarf tomato). The CROP numbers are tomato-
specific and do NOT apply to our 2018 cucumber; the GREENHOUSE EQUIPMENT and
ECONOMICS are generic Dutch-greenhouse and usable to set our twin correctly.

## Validation of our design

- Control = PID: PI heating setpoint (closely followed) + PROPORTIONAL-BAND vents,
  with a dead-zone (vent setpoint = heating setpoint + offset). Exactly our V3.0 design.
- Crop = source–sink with per-organ sink strengths, updated ONCE PER DAY. Exactly our
  Marcelis growth model + daily loop. Independent confirmation.

## Actuator capacities (usable)

- **CO2 dosing capacity `pureCO2cap` = 100 kg/(ha·h)** = 2.78 mg/m²/s. Ours (`psi2`)
  delivers ~1.33 mg/m²/s ⇒ ~2× undersized. Fix: override `psi2 ≈ 27800` in V3.0
  (2.78 × alpha6). Dosing capacity is REDUCED when vents open (table: vent 20%→100%,
  40%→50%, 70%→25% of max).
- **LED lamps**: 2.7 µmol/J; typical 200 µmol/m²/s ⇒ 74 W/m² electric. (matches AiCU ~187.)
- **Ventilation**: ~15% of floor area; lee side first, then wind side.
- **Screens**: TWO — energy screen (scr1, 70% light transmission, closes at night to
  save heat) + blackout (scr2, <1%). Our model has ONE screen (U1). The energy screen
  is a likely contributor to our night heat-loss/undershoot.
- **Slope limits**: heating setpoint ±2 °C/h; CO2 setpoint ±500 ppm/h.

## Control refinements to adopt in V3.0

- **P-band widens** when cold and windy: ~4 °C at 18 °C outside → ~20 °C at 8 °C
  outside (interpolated), + 0.5 °C per m/s of wind above 6 m/s.
- **Vent offset (dead zone)**: ~1 °C at night, 2–3 °C by day.
- **radiationInfluence**: raise heating setpoint with light, e.g. +2 °C as solar goes
  100 → 400 W/m² (light–temperature coupling growers use).
- Setpoints can be time-of-day AND date schedules, and relative to sunrise/sunset.

## Net-profit reward (Stage 3) — cost coefficients

Net Profit [€/m²] = Gains − Fixed costs − Variable costs.

Variable (operational — the agent's trade-offs):
- **Heating**: €0.09 / kWh. `HeatingCosts = Σ(pipe power W/m²)/1000 × 0.09`.
- **Electricity (lamps)**: on-peak (07–23 h) €0.30/kWh, off-peak €0.20/kWh.
- **CO2**: pure CO2 €300/ton = €0.30/kg. `CO2Costs = kgCO2 × 0.30`.

Fixed (per m²·year × fractionOfYear = cropdays/365):
- Greenhouse depreciation/maintenance: **€15.00 / m² / yr**.
- CO2 dosing system: €0.015 / m² / yr per (kg/ha/h) of `pureCO2cap`.
- Lamp depreciation: €0.07 / (µmol m⁻² s⁻¹) / yr × intensity.
- Screen: €1.00 / screen / yr.
- Plant cost (tomato €0.75/plant × plants/m²) — for cucumber use 2018 plant cost.
- Spacing system: €1.50 per spacing change / yr (tomato-specific).

Gains: tomato priced by fresh fruit weight (≥350 g) and DMC (7–8%), × plants/m² — for
CUCUMBER we use the AGC-2018 `Prod_value_cum` / cucumber price from Production.csv.

## Crop-specific (NOT for cucumber, for reference only)

- Dwarf-tomato plant density 56 → 20 plants/m² (potted, spacing over time). Our
  high-wire cucumber is ~2.5 stems/m² — unchanged.
- Quality = DMC + ripeness DVSfruit (≥0.5 sellable).

## Actions

1. DONE-able now: override `psi2 ≈ 27800` in the V3.0 control twin (match 100 kg/ha/h).
2. Investigate night heating undershoot: energy-screen behavior + cover insulation
   (the doc gives no single max heating W/m², so reconcile via the screen + loss terms).
3. Adopt control refinements (P-band widening, day/night vent offset, radiationInfluence,
   slope limits) into `climate_computer.jl`.
4. Stage 3: build the net-profit reward from the coefficients above (cucumber economics).
