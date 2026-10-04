# HANDOFF DOCUMENT — "Greenhouses" project
Prepared 2026-10-04 for migration to a new Claude account.
Source of truth: the project's 19 knowledge docs (notes/01 … notes/19) plus the project settings. Nothing here comes from chat history (see Section 7).

> **FOR OTHER PROJECT CHATS — HOW TO COMPLETE THIS DOCUMENT**
> This draft was written from the 19 project notes only. If you are a chat in the old project with access to earlier conversation history, please add what is missing without deleting existing content:
> 1. Section 2: paste any custom instructions that were ever used (the settings field was empty when this was written).
> 2. Section 6: add personal preferences, phrasing habits, and rules the user stated in chat that are not in the notes.
> 3. Section 7: resolve each flagged gap you can (decisions never written into notes, files uploaded in chat, ensemble-eval result, GitHub repo name, deadlines, paper venue).
> 4. Section 3/4: add decisions or in-progress work that appear in chat but not in the notes.
> Mark every addition with `[ADDED by chat, <date>]` so it can be told apart from the original draft.

---

## 1. Project purpose and intended audience

**Project description (copied verbatim from the project settings, typos included):**
> This project has two main goals.
> 1) To finish a tool for adjusting athe parameters greenhouse+vegetable growth model featuring photosynthesis to make a digital twin.
> 2) To develope an IA agent that controls this greenhouse.

**Plain-language restatement**
- **Goal #1 — Digital twin.** A fast, calibrated simulation of a greenhouse. It couples a Vanthoor-type climate model to FvCB C3 photosynthesis and a cucumber growth model. Its parameters can be inferred from data by name, using t-walk MCMC.
- **Goal #2 — Control agent.** An AI agent that sets the greenhouse climate program (setpoints) to maximize net profit. It is trained on the twin and benchmarked against CEM, MPC, and real growers.

**Audience / context (inferred from the notes, not stated explicitly)**
- Academic research at UNAM (the source folder path is `~/Documents/trabajo/academia/unam/proyectos investigacion/bio/invernaderos`). The user is technically strong: Julia, ODEs, Bayesian inference, RL.
- Outputs are papers and technical reports (`docs/greenhouse_model.{tex,pdf}`, `docs/agent_results.{tex,pdf}`, `docs/greenhouse_control_report.pdf`).
- Calibration data is from the Autonomous Greenhouse Challenge 2018 (AGC-2018, cucumber, six teams). Economics and equipment values come from the AGC-2024 simulator document (dwarf tomato; used only for equipment and cost coefficients).
- Local machine: the desktop device was named "santiagos-laptop-local". The Anthropic account email is antonio.capella.kort@gmail.com.

---

## 2. Complete custom instructions, copied word for word

**The project's custom-instructions field is EMPTY.** The `instructions` value returned by the project settings was an empty string `""`. There is nothing to copy. If you wrote instructions elsewhere (a previous version, another project, or account-level "personal preferences"), they could not be seen. Back them up manually (Section 7).

The working conventions that appear in the notes are summarized in Section 6. They are observed practice, not formal instructions.

---

[ADDED by chat, 2026-10-04] No formal project custom instructions or standing
instructions were given by the user in this chat. The draft's statement that the
field is empty holds for chat history too. What existed were consistent *working
agreements* established through repeated requests (gating long runs, exact run
commands, honest caveats) — these are observed preferences, recorded in Section 6,
not formal instructions. (unsure: whether any account-level "personal preferences"
exist; they were never visible to this chat.)

## 3. Key decisions and outcomes so far

### 3.1 Software architecture and versioning (notes 01–03, 07)
- Original code lives in three Julia folders: `invernaderos2.0` (mature, global-variable-heavy), `invernaderosV2.1` (a refactor that never ran), and `GreenHouse_3.0` (Dash GUI).
- **Decision:** build **V2.2** (`invernaderosV2.2/GreenhouseSim`) from 2.0's proven physics. It has no globals. All 151 constants (115 climate, 36 crop) sit in one concrete `NamedTuple` passed via `p`. Any parameter subset can be inferred by name through `update_params(base, names, x)` (the `merge` pattern).
- **Result:** one `rhs!` call takes 338 ns with 0 allocations. A 2-day `forward_map` takes 5.6 ms and 59 KiB, versus ~19–28 s and 13.33 GiB for the old model. `@code_warntype` is clean.
- **V2.2 is FROZEN** (see its `FROZEN.md`) as the Goal #1 engine. Later realism changes are default-off flags or live in the control layer.
- **V3.0** (`invernaderosV3.0/GreenhouseControl`) is the control/agent layer. It `include`s V2.2 and never copies it, which avoids version drift. A planned monorepo reorg is `GreenhouseCore + Inference + Simulator + Control`.
- The claude.ai project stays "Greenhouses" and spans both goals.

### 3.2 Growth model and calibration (notes 04–06, 12)
- Marcelis-style source–sink partitioning on a daily outer loop. Within-day, the fast ODE runs at a fixed LAI and integrates assimilation `A`. Once a day the growth model spends it on maintenance and growth, sets and harvests fruit cohorts, and updates LAI. This also fixed the 2.0 quirk `LAI = I9`.
- **Calibration** used t-walk MCMC with a multi-objective likelihood on measured climate. It fits weekly yield, leaf/node number, and derived LAI, each with its own inferred sigma. It went through 8 iterations.
  - Beer's-law light interception (`k_ext`) fixed the earlier ~50% overshoot.
  - Topping at `top_day=104` fixed the leaf plateau.
  - Carbohydrate-regulated fruit set produces realistic harvest flushes.
- **Cross-team calibration** (all six AGC-2018 growers, same model and priors) gave yields within ~10%.
  - Development parameters cluster (`node_rate` 0.077–0.099, `set_start_day` 12–16 d), so they transfer across growers.
  - `J_max` (93–198 µmol) and `veg_sink_max` (24–37) vary about 2:1.
  - **Caveat:** that spread is partly confounded by the shared assumption of 2.5 stems/m² and 0.05 m² leaf area. Real per-team densities are still needed.
- **Scientific finding:** under AiCU's CO2 strategy the crop is light-limited, so yield tracks light and `J_max`, not Rubisco capacity (`V_cmax25` is unidentifiable).
- **Cohort / layered canopy (note 12):** leaf-age cohorts, thermal-age fruit cohorts, and per-layer photosynthesis along the light gradient.
  - The layered model fixes more carbon than big-leaf, so the ~15% yield gap was mostly an artifact. It gives 33.5 kg/m², close to the measured ~31, with no recalibration.
  - Optimal density is governed by `k_ext` (3.5 stems at k=0.9; 4.5 at k=0.75). Economics barely move it.
  - Thinning later never pays in the carbon/profit balance.
  - `k_ext` and density are confounded in the data, so density is treated as an agent control variable and `k_ext` should be fixed from the literature (~0.85–0.95).

### 3.3 Physics corrections (notes 07, 09, 10, 11)
- **Leakage:** `nu4` and `nu4CO2` had been inflated ~100× (1e-2 vs Vanthoor's 1e-4) to stabilize the ODE. That caused night undershoot and CO2 that would not enrich.
- **Ventilation bugs:** (1) a spurious `sign(T2 - I5)` factor on the roof-vent rate `f7`; (2) the roof vent `f4` was missing from the sensible-heat vent loss `h7`. The second one caused ~40 °C runaway on calm days. Both are original-model bugs, also present in 2.0.
- With the fixes, leakage was restored to **1e-4**. A venting-day MCMC gave `nu4` ≈ 1.6e-4, with the old 1e-2 excluded.
- **Floor thermal mass:** one dynamic floor state `T5` was added. It is calibrated as a single conductance `k_ground = 4.33` W/m²/K (95% CI [1.45, 7.82]), with `C_soil = 1.5e5` J/m²/K fixed.
- **Heat-balance audit** (`heat_audit.jl`) found FIR-to-sky as the dominant loss (~164 W/m²). Three fixes followed:
  1. Cloud-derived sky temperature (Berdahl-Martin).
  2. A real low-e night screen with `eps_screen = 0.4` and a humidity crack.
  3. VPD-based humidity control (`VPDmin` instead of RH).
- **Result:** heating went from 208 to ~100–129 kWh/m², against a reference of ~1.18 kWh/m²/d.
- Net reward accounting was validated against the real AGC-2018 Reference grower: CO2 0.088 vs 0.079 kg/m²/d and lamp electricity 1.18 vs 1.26 kWh/m²/d.

### 3.4 Control architecture (notes 07, 08, 11)
- **Hierarchical:** the agent acts in setpoint space and a fixed PID layer executes. This was independently confirmed by the AGC-2024 simulator doc.
  - PI heating on air temperature.
  - Proportional-band venting with a dead-zone offset.
  - PI CO2.
  - Rule-based screens and lamps.
- **Action space** grew from 9 to 10 to 11 daily setpoints (all normalized). The set is Tday, Tnight, CO2, light_hours, light_start, VentpBand, ofset, ToutMax, VPDmin, ScreenRad, plus lamp dimming and flue-gas CO2 (inc1). Density is an episode-level choice.
- **Constraints:** slope limits (heat ±2 °C/h, CO2 ±500 ppm/h), temperature ~12–30 °C, CO2 400–1200 ppm, photoperiod 0–20 h.

### 3.5 RL arc (notes 10, 13–19)

| Stage | Result (mean TEST profit, EUR/m2) |
|---|---|
| CEM constant program | −43.8 → −0.29 (big-leaf twin) |
| Cohort twin, CEM | −24.1 → −0.77 |
| Early SAC (24k steps, 3 seeds) | Did not beat CEM (flawed run, see below) |
| SAC long-run, plain replay (15 seeds) | +0.40 ± 1.12; 12/15 beat CEM |
| SAC, demo-anchored replay (15 seeds) | **+1.48 ± 0.68; 15/15 beat CEM** |
| inc1 (lamp dimming + flue-gas CO2), unprotected | −0.13 ± 0.40; 4/15 (broken) |
| inc1b (critic warm-up + decaying BC tether, 100k steps) | **+0.76 ± 0.58; 14/15; deployed seed 104: +1.66 vs CEM 0.12** |
| MPC-distilled SAC (state-conditioned BC from perfect-foresight CEM demos) | **+2.35 ± 0.30; 15/15; deployed seed 202: +2.83** |
| Perfect-foresight ceiling (open-loop CEM per year) | ~3.65–3.70 (loose upper bound) |

**Key lessons**
- Notes 13–14's conclusion that SAC cannot beat CEM was wrong. The causes were training stopped before a late phase transition at ~30–33k steps, model selection on the tuning years, and a corrupt 2023 weather year (all-zero data). The methodology was the finding.
- Fixes that worked:
  - Disjoint TRAIN/VAL/TEST year splits.
  - Validation best-checkpointing.
  - lr 1e-4, target-entropy −3, no random warm-up.
  - Stiff-solver fallback.
  - Demo-anchored replay.
  - Critic-only warm-up of 4000 updates.
  - Decaying BC tether (3.0 → 0 over 40k steps).
- The `forecast_skill` semantics are **0 = PERFECT forecast, 1 = nominal noise, larger = worse**. Earlier notes had this backwards (corrected in note 18).
- The agent **ignores its weather forecast**. The curve is flat: 2.83 / 2.87 / 2.85 / 2.84 / 2.79 across skill 0 to 1.
  - Its gain is a better seasonal and state-conditioned program shape (phenological light budgeting, CO2 ramped with the canopy, temperature at the 16 C cost floor, lamps placed off-peak before 07:00).
  - Decision: **do not extend the forecast horizon `_FC_H` from 3 to 7.**
- **Responsive stomata** (Jarvis, in series with a mesophyll term, flag `stomata_model`) and **diurnal respiration** (`resp_diurnal`) were validated but are second-order in light-limited autumn. Default stays `stomata_model = 0`. Revisit for a high-light spring/summer season by re-running `recalibrate_stomata.jl`.
- **CO2 enrichment is a net loss in autumn.** Pg rises only ~7% by 1200 ppm and profit falls.
- **vs real growers (note 19):** the agent grows less (22.8–28.3 vs 36–51 kg/m2) but spends far less. Its profit is nearly flat versus electricity price, while growers are steeply price-sensitive. The crossover is ~0.10–0.12 EUR/kWh. At cheap 2018-era power (0.05) the best growers (Sonoma 16.4, Reference 14.0, iGrow 13.9) out-earn the agent (9.9). Verdict: **the agent is price-robust, not a better grower.**
- **OOD finding:** on the higher-yield Reference crop the agent goes erratic in the final weeks (night setpoint above day, temperature crash). It is out-of-distribution on the crop-state observation (cumulative yield), not a twin bug. The fix is a short fine-tune on the matched crop.

---

### Additions from chat history

[ADDED by chat, 2026-10-04] **Stomata-modeling fork (decided via an explicit
multiple-choice question).** Chosen: add responsive stomata as an *opt-in flag
inside the frozen V2.2 engine* (default off), using a *Jarvis multiplicative*
conductance. Rejected alternatives, with reasons:
- Fork a new version folder (V2.3) — rejected (version-folder proliferation was
  already a documented pain; one frozen engine is better).
- Edit V2.2 in place without a flag — rejected (would break bit-for-bit
  reproducibility of all prior results).
- Ball–Berry / Medlyn (A-coupled) conductance — rejected because gs depends on
  assimilation, needing an inner fixed-point solve on every call in the ODE hot
  loop (cost + convergence-stability risk on stiff days).

[ADDED by chat, 2026-10-04] **Series-mesophyll refactor of the stomata model.**
The first version let the Jarvis term set the *total* conductance; a
calibration-neutral fit then pushed `gs_max` to ~0.57 (unphysically high for a
total conductance). Decision: make Jarvis set only the *stomatal* component in
series with a fixed mesophyll conductance `gs_mesophyll = 0.25` (mirrors the frozen
`g_eff = 0.3*0.25/(0.3+0.25)`), which caps total conductance physically.
Recalibrated `gs_max = 1.2469` (calibration-neutral: matches the frozen twin's
seasonal gross assimilation at ambient CO2). Full data-refit of stomata via t-walk
was rejected for now (no leaf gas-exchange data). Validation: `validate_stomata.jl`
TEST 1-4 all PASS (bit-for-bit off; monotone; anchor ratio 1.00).

[ADDED by chat, 2026-10-04] **"Bigger campaign vs longer runs vs architecture?"**
Resolved with evidence from the per-seed VAL curves: *longer runs help* (good
seeds' VAL still climbing at the 60k cutoff; the 11-action problem's phase
transition moved to ~50k), *architecture does not* (good seeds reach good policies,
so it's a training-dynamics problem not a capacity one), *more seeds is secondary*.

[ADDED by chat, 2026-10-04] **Direction fork after the "agent ignores its forecast"
finding (decided via an explicit multiple-choice question).** Chosen: "MPC teaches
RL (both)" — build the perfect-foresight optimizer, then distill its
forecast-driven schedules into SAC as state-conditioned demonstrations. Deferred
(not rejected) alternatives: (A) a receding-horizon MPC controller on its own; (B)
fix RL forecast use directly via a longer horizon + forecast-aware training; (C)
conclude the RL line and only write it up.

[ADDED by chat, 2026-10-04] **MPC->RL distillation outcome (the chat's main new
result).** Deployed MPC-distilled agent = +2.83 EUR/m2 (15/15 seeds beat CEM; mean
+2.35 +/- 0.30), ~78% of the ~3.65 ceiling, vs reactive +1.66. Mechanism: expert
demos from the open-loop CEM optimizer per TRAIN year (warm-started from SAC),
imitated with a *state-conditioned* behavior-cloning term that replaced the
forecast-blind constant-program tether. Smoke test alone (2 demo years) already
gave TEST 1.96 from the behavior-cloning warm-start, which also makes the
imitation performance a *floor* for the campaign.

[ADDED by chat, 2026-10-04] **Report decision.** Rather than update the stale
`docs/agent_results.{tex,pdf}`, a *new consolidated control-agent report* was
written and delivered: `docs/greenhouse_control_report.{tex,pdf}` (5 pp). It
covers the full ladder, the forecast-flat finding, the AGC comparison, the
electricity-price sweep, and (added after reviewing the plots) an
"Out-of-distribution diagnostic" paragraph. The To-Do is led — at the user's
explicit request — by price-adaptive agents.

[ADDED by chat, 2026-10-04] **Reorg architecture recommendation (accepted in
principle; execution pending GitHub credentials).** One git monorepo with a shared
`GreenhouseCore` engine + `Inference` + `Simulator` + `Control` layers; version
milestones by *git tag/branch, not version folders*; per-goal chats are fine as
long as they share the one git repo; a *private* GitHub repo recommended; sequence
= reorg first, then Goal 1 (inference + probabilistic forecast), with Goal 2
(simulator) in parallel/after.

## 4. Work in progress and next steps

**Status:** the RL line is **paused** at the consolidated report (`docs/greenhouse_control_report.pdf`, delivered). Goal #1 is considered effectively complete (notes 06, 09), with a probabilistic production forecast still to do.

**Next steps, in the order written in note 19**
1. **Git monorepo reorganization** (GreenhouseCore + Inference + Simulator + Control). **Blocked: waiting for the user's GitHub credentials.**
2. **Goal #1 follow-through:** run t-walk parameter inference and build a **probabilistic 3–5 week production forecast**.
3. **Price-adaptive agents** (the lead To-Do in the report). Put electricity and energy prices in the observation and train across price regimes. A retrain at cheap power should make the agent light up.
4. **Fine-tune on the matched crop** to fix the OOD behavior on the higher-yield Reference crop.

**Other open items from the notes**
- Real per-team plant densities (AGC protocol docs). Re-run the cross-team calibration with them.
- Fix `k_ext` from the literature (~0.85–0.95).
- Optional humidity/disease penalty so thinning can pay for non-carbon reasons.
- Update `docs/agent_results.{tex,pdf}` with the demo-anchored and ensemble numbers.
- Ensemble evaluation (`ensemble_eval.jl`, top-5 action-average) was pending in note 15. Its outcome is not recorded.
- Update tariffs. The AGC-2018/2024 prices are pre-energy-crisis and probably low.
- Actuator-freedom candidates (ranked for autumn): second screen, continuous dimming (done), minimum pipe temperature, flue-gas CO2 (done), richer temperature schedule, CHP + buffer + electricity selling.
- Only the autumn crop (planted Aug 14, the lowest-light slot) is covered. Three NL planting seasons exist. A spring/summer season would need stomata recalibration and a re-baseline of CEM and RL.
- Check whether VAL fully plateaued at 100k steps. Seed 112 still lost (−0.11).
- The big-leaf twin's absolute yield (19 kg/m2) is unverified against a matched season. The AGC-2018 weather is needed in `load_weather` form.
- Minor latent bug: the `env_step!` success return omits `failed` while the failure return includes it. Fold `failed = false` into the success return.
- V2.2 `assimilation` still takes `LAI = I9` in the non-growth path (a faithful 2.0 quirk). The growth path uses the real LAI.
- `test/benchmark.jl` needs re-baselining (its reference assumed `nu4` = 1e-2).
- Possible speedup: cached-index or gridded interpolant, since interpolation dominates `rhs!` (~316 of 338 ns).

---

### Additions from chat history

[ADDED by chat, 2026-10-04] At stop, the project is **mid-migration to a new Claude
account** (this handoff is the current activity). No compute was running when we
stopped.

[ADDED by chat, 2026-10-04] **Done in this chat but still leaving open follow-ups:**
- Plots generated (`agent_vs_grower_native.png`, `agent_vs_grower_matched.png`);
  the OOD anomaly was *confirmed* (native run smooth, matched run erratic). The
  recommended **fine-tune on the matched Reference crop** to fix it is NOT done.
- The electricity-price sweep was a *re-scoring only* (agent policy is fixed); a
  **retrain at cheap power (~0.06 EUR/kWh)** for a fair cheap-power comparison is
  NOT done.
- `docs/agent_results.{tex,pdf}` is still **not updated** and is now stale; the
  new consolidated report supersedes it but the old file remains.

[CONFLICT] Draft Section 4 lists "Update docs/agent_results with the demo-anchored
and ensemble numbers" as an open item, implying that is the reporting path. Chat:
that doc was NOT updated; instead a new consolidated report
(`docs/greenhouse_control_report.pdf`) was written and delivered. Both are true —
the old item is still open and the new report now carries the numbers — but the
current report is the consolidated one, not an updated `agent_results`.

## 5. Knowledge base — every uploaded file

**No uploaded files and no sync sources exist.** The project file list was empty (`files: []`). The knowledge base is **19 Markdown docs** (~29 KB total). Dates are creation dates.

| Path | One-line description |
|---|---|
| notes/01-model-architecture.md (2026-09-08) | Map of the Julia code (2.0 / V2.1 / 3.0 GUI): shared Vanthoor+FvCB model, state vector, how to run each, and the missing growth model. |
| notes/02-performance-assessment-v2.1.md | Diagnosis that V2.1 is not faster because the RHS still reads non-const globals; design for a type-stable NamedTuple + merge inference. |
| notes/03-v2.2-fast-forward-map.md | Build log and benchmarks of V2.2 (338 ns RHS, 0 alloc, 5.6 ms 2-day forward map), file layout, run commands, bring-up gotchas. |
| notes/04-cucumber-growth-model.md | Design of the Marcelis source–sink cucumber growth model, its daily loop, parameters, and coupling into V2.2. |
| notes/05-growth-calibration.md | Eight-iteration t-walk calibration of growth parameters on AiCU yield, leaf, and LAI; identifiability findings. |
| notes/06-cross-team-calibration.md | Calibration across all six AGC-2018 growers, the parameter table, and the density-assumption confound. |
| notes/07-control-agent-architecture.md | V2.2 frozen / V3.0 control layer decision, setpoint+PID architecture, Stage-1 validation, and the leakage-inflation root cause. |
| notes/08-agc2024-reference-values.md | Equipment capacities, control refinements and the cost/net-profit coefficients from the AGC-2024 simulator doc. |
| notes/09-ventilation-fix-and-leakage-validation.md | Two original-model ventilation bugs, leakage re-validation at 1e-4, and the calibrated floor thermal-mass state. |
| notes/10-stage3-rl-design-and-economics.md | RL MDP design, algorithm landscape, net-profit reward, heat-balance audit/FIR fix, weather and forecast design. |
| notes/11-control-surface-and-realism-fixes.md | Full action space, density schedule, the three control realism fixes, and baseline economics. |
| notes/12-vertical-canopy-and-density.md | Cohort/layered canopy model, density-vs-k_ext optimum, thinning conclusion, and the k_ext/density confound. |
| notes/13-rl-env-foundation-and-cem-baseline.md | VPD control, weather episodes and noisy forecast, stepping env, CEM baseline, first SAC agent. (Later conclusions superseded by note 15.) |
| notes/14-cohort-twin-and-rl-vs-cem.md | Cohort crop in the env, lamp-blindness fix, multi-seed result that SAC did not beat CEM. (Superseded by note 15.) |
| notes/15-rl-long-run-and-data-fix.md | Corrected methodology (clean data, splits, best-checkpoint, demo-anchored replay): SAC beats CEM 15/15; 2023 bad-data fix. |
| notes/16-responsive-stomata-and-diurnal-respiration.md | Opt-in Jarvis stomata and diurnal respiration; validation, recalibration, and why they are second-order in autumn. |
| notes/17-inc1-campaign-and-warmstart-fix.md | Lamp dimming and flue-gas CO2 increment; the instability diagnosis and warm-start fix; 14/15 seeds beat CEM. |
| notes/18-mpc-ceiling-and-forecast-value.md | Perfect-foresight MPC ceiling (3.70), correct forecast_skill semantics, and the finding that the agent ignores its forecast. |
| notes/19-mpc-to-rl-distillation.md | MPC→RL distillation (2.83), the flat forecast curve, comparison to real growers, price sweep, OOD diagnostic, and next steps. |

**Reading order for a newcomer:** 01 → 03 → 04 → 06 → 09 → 07 → 10–11 → 15 → 17 → 18 → 19 (skim 13–14 only for history).
**Note numbering:** there is no `notes/00` and nothing missing in the sequence. Notes 07 and 08 were created in reverse order (08 first), so note 07 cites note 08.

**Files that exist only on the user's computer or in the repo (not in the project):** the Julia source folders, `configfiles/*.json` (including `climate_parameters.json`, `constants_cropphoto.json`, `cucumber_growth.json`), AGC-2018 CSVs (`Production.csv`, `CropManagement.csv`, `Greenhouse_climate.csv`, `meteo.csv`), `dataset_meteo_holanda.csv`, trained actors (`test/sac_actor_best.jls`, `test/sac_actor_mpc_best.jls`, `test/campaign_actors/`), `docs/*.tex|pdf`, and the AGC-2024 "Part B" PDF. Back them up separately.

**Important scripts referenced:**
- V2.2: `test/benchmark.jl`, `test/run_inference.jl`, `test/run_inference_venting.jl`, `test/run_growth.jl`, `test/check_growth.jl`, `test/calibrate_growth.jl`, `test/calibrate_all_teams.jl`, `test/validate_stomata.jl`, `test/recalibrate_stomata.jl`.
- V3.0:
  - Control and physics: `test/run_control.jl`, `test/heat_audit.jl`, `test/sweep_leakage.jl`, `test/validate_leakage_agc.jl`.
  - Optimization and RL: `optimize_cem.jl`, `test/train_sac.jl`, `test/train_sac_long.jl`, `test/train_sac_mpc.jl`, `test/train_sac_multiseed.jl`.
  - Evaluation: `test/ensemble_eval.jl`, `test/mpc_ceiling.jl`, `test/forecast_value.jl`, `test/gen_mpc_demos.jl`, `test/compare_vs_agc.jl`, `test/plot_actions.jl`, `test/diagnose_policy.jl`, `test/diagnose_limitation.jl`.
- Source files: `src/climate.jl`, `photosynthesis.jl`, `growth_cohort.jl`, `forward_map.jl`, `parameters.jl`, `weather.jl`, `weather_episodes.jl`, `rl_env.jl`, `climate_computer.jl`, `economics.jl`, `twin.jl`.

---

### Additions from chat history

[ADDED by chat, 2026-10-04] **No data files were uploaded by the user in this chat.**
The user pasted *terminal outputs* (campaign summaries, test tables, the
forecast/price-sweep results, plot observations) and, in the final message, the
text of the handoff document itself. New files produced in this chat that are not
in the draft's list:
- `docs/greenhouse_control_report.tex` / `.pdf` — the consolidated control-agent
  report (delivered to the user and committed to `docs/`).
- `test/mpc_demos.jls`, `test/mpc_demos_smoke.jls` — forecast-driven expert
  demonstrations (1800 / 200 transitions) for MPC->RL distillation.
- `test/sac_actor_mpc_best.jls` (in draft) plus `test/campaign_actors_mpc/`
  (per-seed MPC-distilled actors) and `test/campaign_manifest_mpc.csv`.
- `test/figures/agent_vs_grower_native.png`, `agent_vs_grower_matched.png`.
- Run logs: `sac_long_inc1.log`, `sac_long_inc1b.log`, `sac_mpc.log`,
  `gen_demos.log`, `mpc_ceiling.log`.
- Config edits (Mac): `constants_cropphoto.json` gained `stomata_model`=0,
  `gs_max`=1.2469, `gs_mesophyll`=0.25, `gs_I_half`, `gs_D0`, `gs_Ca_ref`,
  `gs_Ca_scale`, `gs_min`; `cucumber_growth.json` gained `resp_diurnal`=0.
  `.bak_prestomata` / `.bak_preseries` backups were made.

## 6. Recurring context, preferences, terminology, and style rules

*Everything below is inferred from the notes. Chat history could not be seen, so these are conventions of the work, not confirmed personal preferences.*

**Terminology**
- **Twin** is the simulator (V2.2 engine + V3.0 control). **Forward map** is parameters → simulated observations. **Frozen** means V2.2's validated engine, not to be edited except through default-off flags.
- **CEM** is the Cross-Entropy Method, and "CEM constant" is the tuned constant-season program baseline. **SAC** is Soft Actor-Critic. **MPC ceiling** is the perfect-foresight upper bound. **Demo-anchored replay**, **warm-start**, **BC tether**, and **distillation** are the RL stabilization tricks.
- **TRAIN/VAL/TEST** are disjoint weather-year splits. Weather bank: 30 usable years (31 minus the corrupt 2023).
- **Autumn crop** is planted Aug 14, light-limited. **AiCU, Croperators, DeepGreens, Reference, Sonoma, iGrow** are the AGC-2018 teams.
- Parameter names:
  - Growth: `J_max`, `V_cmax25`, `k_ext`, `asrq`, `node_rate`, `veg_sink_max`, `set_rate`, `set_start_day`, `LAI_max`, `Wf_max`, `top_day`, `stem_density`.
  - Climate: `nu4`, `nu4CO2`, `k_ground`, `C_soil`, `gamma4`, `psi2`, `eps_screen`.
- Controls **U1…U12**: thermal screen, pad-fan, mech cooling, heater, shade screen, side vents, forced vent, roof vents, fog, CO2, pipe setpoint, lights.
- Units: profit in **EUR/m2** (per season, ~100 days). Yield in **kg FW/m2** (dry mass / DMC 0.04). Price 0.889 EUR/kg FW. Electricity 0.30 peak (07–23 h) and 0.20 off-peak.
- Heating 0.09 EUR/kWh, CO2 0.30 EUR/kg.
- Reward is scaled x10 for SAC learning while reported euros stay true.

**Working preferences visible in the notes**
- Rigor over headline wins. Multi-seed replication, disjoint held-out evaluation, reward audits, and explicit notes retracting earlier wrong conclusions (13–14 → 15) are standard practice.
- Validation must be explicit (PASS/FAIL tests, regression against 2.0, heat-balance audits within ~1%). Flags default OFF to keep the frozen engine bit-for-bit reproducible.
- Honest framing: say when results are not apples-to-apples, and state caveats (density confound, OOD, price-robust ≠ better grower).
- Fast, type-stable Julia only on the hot path. Ordinary Julia is fine for the daily growth step.
- Julia-native tools preferred (Flux for SAC, kept out of the core package). The user had no preference between Julia RL and a Python Gym wrapper.
- Notes are numbered Markdown files in `notes/`, structured as headline result → method → table → caveats → next. Reports are LaTeX/PDF in `docs/`.
- Environment habits: always `--project=.`, reuse 2.0's `Project.toml`/`Manifest.toml`, and use `OrdinaryDiffEq` (not `DifferentialEquations`). Use `Tsit5()` (stiff solvers TRBDF2/Rodas as fallback). The DataInterpolations API changed: use `extrapolation=ExtrapolationType.Constant`.
- Work on the user's files on their computer. **`device_commit_files` sometimes reports success but does not write**, so always verify on disk with grep and retry or edit in place.

---

### Additions from chat history

[ADDED by chat, 2026-10-04] **Hard execution constraint (confirmed in chat):** the
assistant cannot run Julia — neither the cloud container nor the mounted device VM
(`device_bash`) has Julia. The **user runs ALL Julia on the real Mac and pastes the
output**. So: prepare scripts and give exact run commands; never assume you can
execute them. (The device shows as `santiagos-laptop-local` via get_device_info but
the Mac prompt is `ACK@MacBook-Air-4`, home `/Users/ACK/...`.) Note: the starter
prompt (Section 8) already encodes the related rules "ask me before starting any
long compute run" and "I'll work on my own computer's files" — this confirms them.

[ADDED by chat, 2026-10-04] **Smoke-test-then-launch pattern.** For every long run,
the user ran a short smoke test first and inspected it before committing hours of
compute. Always provide a smoke mode and wait for confirmation before the full run.

[ADDED by chat, 2026-10-04] **Confirm before long compute; execute in order;
mechanism first.** Representative user instructions (verbatim):
- "respond thsi before moving to the RL." (gate a step on an answer first)
- "Lets leave the extintion and plant density like it is, and move to the RL agent."
- "do these thing in order."
- "We will do all 3, but first confirm the mechanism."
- "It is clear that the control problem is harder."

[ADDED by chat, 2026-10-04] **Likes sensitivity analyses and visual diagnostics.**
- "can you re-run the test with the cheaper electricity prices and see what happens?"
- "I would also like to look at the plots of some variables to see how different
  the actions were."
- "Make the report and mention in the ToDo the ideas of agents that respond to
  electricty, gast, etc prices." (origin of the price-adaptive-agents To-Do — a
  *user-requested* direction.)
- "add the comment, I wil give you the github credentials later on."

[ADDED by chat, 2026-10-04] **Honest, non-overclaiming framing is valued.** The user
engaged positively with caveated conclusions and asked for the test that could
overturn a flattering result (the cheap-power sweep that showed the agent is
"price-robust, not a better grower").

[ADDED by chat, 2026-10-04] **Language note (factual, possibly a style cue):** the
user is at UNAM (Mexico) and writes English with Spanish-influenced phrasing and
typos — e.g., "IA agent" = AI agent (Spanish word order: *agente de IA*),
"develope", "extintion", "preasure", "gast". Read intent generously through typos.

[ADDED by chat, 2026-10-04] **Terminology added in this chat:** "inc1" = increment 1
(continuous lamp dimming + flue-gas CO2); "inc1b" = the warm-start-fixed inc1
campaign; "MPC->RL distillation"; "expert demos" / "state-conditioned behavior
cloning"; "matched crop" vs "native crop" (the AGC comparison runs the agent on a
real team's calibrated crop vs its own trained crop); "price-adaptive agents"
(prices in the observation, trained across price regimes).

## 7. Anything that could not be recalled or accessed — BACK THESE UP MANUALLY

1. **Custom instructions:** the field is **empty**. If instructions existed anywhere else, copy them manually.
2. **Chat history:** earlier conversations in this project could not be seen, only the 19 saved notes. Reasoning, rejected ideas, and user quotes that never made it into notes are lost unless chats are exported.
3. **Project memory:** the memory tool returned "not available in this session." Open claude.ai → Project → Memory and copy any entries manually.
4. **Uploaded files and PDFs:** none were present. If PDFs (for example the AGC-2024 Part B PDF) or datasets were uploaded in the past, they are not in the file list. Re-upload from the originals.
5. **Everything on the user's computer:** source code, JSON configs, datasets, trained actors, LaTeX/PDF reports (see Section 5). The notes name them but do not contain them.
6. **Numbers that could not be verified:** the notes disagree or are incomplete in a few places.
   - Note 13 reports CEM at −0.29 on the big-leaf twin and note 14 at −0.77/−0.80 on the cohort twin. These are different twins. Later campaigns use a CEM baseline of −0.37 (10-action) and 0.12 (inc1, density 2.40).
   - The ensemble evaluation result (note 15) was not recorded.
   - The MPC ceiling appears as a range (3.60–3.70).
   - GitHub repo name and credentials are not in the notes.
7. **Account-level settings:** name, connectors, installed plugins, skills, and scheduled tasks were not reviewed. Whether any scheduled tasks exist was not checked.
8. **Unstated assumptions that could not be confirmed:** the user's role and institution beyond the UNAM folder path, the paper venue, and deadlines. None are in the notes.

---

### Additions from chat history (gap resolutions)

[ADDED by chat, 2026-10-04] **Ensemble-eval result — RESOLVED (from the
pre-compaction summary in this chat):** "Ensemble(5): mean +2.75 vs best-single
+2.89, better worst year (+1.22 vs +1.08)." Reading: a top-5 action-averaged
ensemble did **not** beat the best single seed on the mean, but had a better
worst-year (more robust). Context: the 10-action demo-anchored campaign era
(note 15), before inc1. (unsure: whether +2.75/+2.89 are absolute profit or
advantage-over-CEM — the summary does not say; magnitudes don't match the note-15
"+1.48" advantage, so treat the exact basis as unconfirmed.)

[ADDED by chat, 2026-10-04] **MPC ceiling range — RESOLVED:** two independent runs
gave mean 3.70 and 3.60 (CEM stochasticity), so ~3.65 is a stable estimate, not a
disagreement.

[ADDED by chat, 2026-10-04] **GitHub repo name/credentials — still UNRESOLVED.**
The user explicitly deferred: "I wil give you the github credentials later on." My
recommendation in chat was a *private* repo named `greenhouses`; not confirmed.

[ADDED by chat, 2026-10-04] **Deadlines and paper venue — still UNRESOLVED.**
Neither was mentioned anywhere in this chat.

[ADDED by chat, 2026-10-04] **Project memory — not checked in this chat either**
(unsure; the memory tool was not used here).

## 8. Ready-to-paste starter prompt for the first chat in the new project

(Also saved separately as `starter_prompt_new_project.md` in this folder.)

> I'm continuing a research project called "Greenhouses". I've uploaded a handoff document and the project notes (notes/01–19). Please read the handoff first, then notes 03, 09, 15, 17, 18 and 19, before answering anything.
>
> **Summary:** Goal 1 is a fast, calibrated Julia digital twin of a cucumber greenhouse (Vanthoor climate + FvCB photosynthesis + Marcelis source–sink growth, V2.2 frozen, calibrated on AGC-2018 for six teams, with t-walk inference by parameter name). Goal 2 is an RL control agent that sets climate setpoints to maximize net profit (V3.0 control layer). The best agent is an MPC-distilled SAC at +2.83 EUR/m2 on held-out weather, against +0.12 for the tuned constant CEM program and a perfect-foresight ceiling of ~3.65. It ignores its weather forecast, is price-robust but not a better grower than the best AGC growers at cheap power, and is out-of-distribution on the higher-yield Reference crop. The RL line is paused.
>
> **Please:** (1) tell me in your own words what you understand the state of the project to be and flag anything in the notes that looks inconsistent or missing; (2) propose a concrete plan for the next three steps, in this order: (a) the probabilistic 3–5 week production forecast with t-walk inference (Goal 1), (b) price-adaptive agents trained across electricity-price regimes, (c) a fine-tune on the matched Reference crop to fix the out-of-distribution behavior. Note that the git monorepo reorganization is blocked until I provide GitHub credentials.
>
> **Ground rules:** keep the V2.2 engine frozen, with new physics as default-off flags; use disjoint TRAIN/VAL/TEST weather years and multi-seed evaluation; report caveats and retractions explicitly; write results as numbered Markdown notes in the project style (headline → method → table → caveats → next); ask me before starting any long compute run. I'll work on my own computer's files, so tell me exactly what to run when you can't.
>
> Remember `forecast_skill`: 0 = perfect forecast, higher = noisier.
