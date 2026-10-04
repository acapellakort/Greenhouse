# 16 — Twin realism increment: responsive stomata + diurnal respiration

Goal: make the twin's response to the agent's humidity/CO2/temperature levers
physically correct, so the RL reward signal is trustworthy. Two additions, both
**opt-in flags that default OFF**, so the frozen V2.2 engine stays bit-for-bit
reproducible. No new version folder — the frozen engine is extended in place
with default-off switches (clean A/B for the paper).

**CONCLUSION (validated): in the light-limited autumn scenario responsive
stomata are correct-but-second-order. Keep `stomata_model = 0` as the autumn
default. The code is validated and physical; it will matter in a high-light
spring/summer season, where it should be recalibrated for that season.**

## What changed (V2.2 GreenhouseSim + V3.0 control layer)

**Stomata — Jarvis stomatal conductance in series with a fixed mesophyll term.**
`photosynthesis.jl` gained `jarvis_gs` (the stomatal component) and a
`stomatal_conductance` dispatcher. `assimilation` takes a `VPD` keyword and calls
the dispatcher instead of the hardcoded `g_eff()`. With `stomata_model = 0`
(default) it returns constant `g_eff` and ignores VPD — identical to frozen.
With `= 1`:

    gs_stom = max(gs_min, gs_max · fI(light) · fD(VPD) · fC(CO2))
    fI = Iw/(Iw+gs_I_half); fD = 1/(1+VPD/gs_D0); fC = 1/(1+max(ppm-gs_Ca_ref,0)/gs_Ca_scale)
    total g = (gs_stom · gs_mesophyll) / (gs_stom + gs_mesophyll)     # series

The **series structure** mirrors the frozen `g_eff = 0.3·0.25/(0.3+0.25)` and
caps total conductance at the mesophyll value (~0.25) no matter how wide the
stomata open — so `gs_max` stays a physical stomatal maximum rather than a fitted
total. Temperature response is left to FvCB (Vcmax/Jmax), not double-counted.
Calibrated params: `gs_mesophyll 0.25, gs_max 1.2469` (calibration-neutral, see
below), `gs_I_half 50, gs_D0 1.0, gs_Ca_ref 400, gs_Ca_scale 900, gs_min 0.01`.

**VPD threading.** Climate RHS passes its VPD straight into the drawdown
`assimilation` call; `_pg_matrix` (rl_env.jl) carries a 4th row = leaf-air VPD
[kPa] into `daily_assimilation_layered`. Both the ODE CO2 balance and the layered
growth see the same responsive conductance.

**Diurnal respiration.** `grow_cohort!` takes an optional `q10_resp` override;
the env passes the Jensen-correct mean of `Q10^((T_hour−25)/10)` over the day's
hourly temperatures when `resp_diurnal = 1` (default 0). Fixes the convexity bias
that underestimates warm-night maintenance respiration. The existing respiration
model was already sound (Q10 × per-organ biomass, growth resp via asrq) — this is
a refinement.

Flags are Float64 entries (type-stable NamedTuple): `stomata_model`,
`gs_*`, `gs_mesophyll` in constants_cropphoto.json; `resp_diurnal` in
cucumber_growth.json. Toggle at runtime via `update_params`.

## Validation (test/validate_stomata.jl — all PASS, series version)

- TEST 1 flag=0 bit-for-bit: max |A(VPD)−A(noVPD)| = 0; gs == g_eff.
- TEST 2 monotone (TOTAL conductance): VPD 0.2→3.0 kPa gs 0.161→0.088; light
  0→600 gs 0.010→0.149 (saturates at the mesophyll cap, as intended); CO2
  400→1400 gs 0.136→0.091.
- TEST 3 anchor: gs(ref) 0.1364 vs g_eff 0.1364, ratio 1.00.
- TEST 4 end-to-end: 0 solver fails both flags.

## Recalibration (calibration-neutral, 1 DOF) — done

`test/recalibrate_stomata.jl` tunes ONLY `gs_max` so the responsive twin matches
the frozen twin's seasonal gross assimilation at **ambient CO2 (400 ppm)** (where
fC=1, pinning intrinsic productivity). Result: `gs_max = 1.2469`. NOTE the
bracket margin was thin — `gs_mesophyll = 0.25` forces the stomatal component
near its cap to match frozen productivity; a slightly higher `gm` would sit more
comfortably. That's a spring/summer refinement.

After recalibration the isolated enrichment/VPD realism penalty on the CEM
program is **−0.04 EUR/m2** — negligible. The earlier "−0.5 to −0.7" penalty
(TEST 4 at the un-recalibrated gs_max=0.33) was almost entirely a too-low
season-mean conductance, NOT CO2 closure.

## Why it washes out: autumn is light-limited (test/diagnose_limitation.jl)

Part A (frozen model), season Pg & profit vs CO2 setpoint:

    CO2 ppm    Pg g/m2   profit EUR/m2
      400      1858.0      -0.52
      700      1862.9      -0.54
     1000      1949.6      -2.26
     1200      1989.2      -2.16

Pg barely moves 400→700 (+0.3%), only +7% by 1200 ppm, and **profit gets worse**
as CO2 rises (cost >> assimilate gain). So (a) the crop is light-limited →
conductance/CO2 have little leverage on assimilation → responsive stomata are
second-order; and (b) **CO2 enrichment is a net loss in autumn**, independently
confirming why the CEM optimum kept CO2 modest. (Part B leaf-level split is
formatted with unscaled mol/m2/s values ~0; the CO2 sweep carries the argument.)

## Downstream

For autumn: `stomata_model = 0` stays the default; the RL work continues on the
frozen (constant-g) twin — no CEM/campaign invalidation. For a future high-light
season: flip the flag, re-run `recalibrate_stomata.jl` for that season's weather
(and consider raising `gs_mesophyll`), then re-baseline CEM + RL on the
responsive twin.

## Minor latent issue

`env_step!` success return omits `failed` while the failure return includes it
(mismatched info fields). Harmless today; fold `failed = false` into the success
return next time rl_env.jl is edited.
