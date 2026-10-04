# 11 — Control surface (agent action space) and control-layer realism fixes

The RL agent programs SETPOINTS (like a grower's climate computer); low-level PI/
proportional controllers turn them into the 12 actuators U1..U12 driving the V2.2
physics. This note records the full control surface and the realism fixes made while
building it. All of this lives in the control layer (V3.0 climate_computer.jl / twin.jl);
climate_rhs stays byte-identical except the calibration-safe eps_screen addition.

## Action space (chosen: Extended-7 + humidity + screen + density)

Daily setpoints (per day):
- `Tset_day`, `Tset_night` [C] — PI heating (pipe temp) and, via the vent band, cooling.
- `CO2_set` [ppm] — PI CO2 dosing (valve U10).
- photoperiod / `light_hours` — lamps U12 (dominant electricity cost).
- `VentpBand` [C] — proportional vent band (roof+side).
- `ofset` [C] — vent dead-band above the heating setpoint.
- `ToutMax` [C] — energy screen also closes below this outside temp.
- `ScreenRad` [W/m2] — radiation below which the energy screen closes at night.
- `RHmax` [%] — humidity setpoint: above it the controller dehumidifies (vent + crack screen).

Episode-level (occasional, not daily):
- **plant density schedule** — set at planting, changed once or twice (interplanting).
  `density_schedule([(0,2.0),(30,2.75),(60,3.0)])`; scales veg_sink_max, set_rate, LAI_max
  in grow!; plant cost scales with density. Effect: higher density -> higher LAI -> more
  yield but DIMINISHING (canopy closes ~LAI 2.5); dense-from-start beats late interplanting.

Couplings the agent must learn: CO2xventing (dosing while venting wastes CO2);
lighting x heating (lamps add ~29 W/m2 heat but cost electricity — pure waste in summer);
temperature x yield x cost.

Constraints to enforce in the env: slope limits (heat +-2 C/h, CO2 +-500 ppm/h), physical
bounds (temp ~12-30 C, CO2 ~400-1200 ppm, photoperiod 0-20 h), vent-setpoint = heat + ofset.

## Control-layer realism fixes (all done, audited)

A season air heat-balance audit (V3.0 test/heat_audit.jl, consistency ~1%) drove three fixes:

1. **Cloud sky** (`load_weather(...; sky_from_clouds=true)`): the twin radiated to a fixed
   clear -0.4 C sky; cloud/RH-derived sky temperature (Berdahl-Martin) cut FIR 394->230 and
   heating 208->140 kWh/m2.
2. **Energy screen** (climate_computer.jl `screen_control`): a real night-closing screen with
   a humidity crack; and the screen-sky FIR term was full-emissivity (epsil5=1, a shade screen,
   RADIATING MORE when closed). Added `eps_screen=0.4` (aluminized low-e) to r12 only
   (calibration-safe, U1=0 on the venting-calibration day). FIR r12 152->61, heating 140->100.
3. **Humidity control** (`humidity_vent`): RH now actively bounded (vent + crack screen when
   RH > RHmax) instead of floating.

Floor thermal mass was ruled OUT as a heating driver (audit: floor on/off 212 vs 214).

## Baseline economics (fixed setpoints, 1998 100-day season, after fixes)

Net -33 EUR/m2: revenue ~17 (18.7 kg @ 0.889), heating ~9 (100 kWh), CO2 ~3.2, electricity
~33 (118 kWh lamps), fixed ~4.7. The loss is dominated by heavy summer lighting and modest
yield -- the levers the agent will pull. Heating now matches reference physics.

## Next: the RL environment (Stage 3b)

Stepping env around simulate_twin: reset! (sample a complete-season weather year), step!
(one day: decode setpoint action -> simulate day -> grow -> daily-profit reward). Obs =
crop + climate + noisy 3-day forecast (note 10). Density = per-episode action. Then a
gradient-free baseline optimizer, then SAC.
