# GreenhouseSim V2.2 — fast, type-stable forward map

A rewrite of the greenhouse **climate + FvCB photosynthesis** forward map,
built from the proven 2.0 physics, engineered for **speed** (for inference and
agent training) while keeping the 2.0 feature of **inferring any parameter by
name** — but *without* global variables.

## Why it's fast (what changed vs 2.0 / V2.1)

2.0 and V2.1 defined every physical constant as a **non-`const` global**
(`eval(:(global …))`). In Julia that makes the compiler treat each one as
`Any`: boxing, dynamic dispatch, no optimization — the cause of the ~13 GiB of
allocations in the old README. V2.2 removes all of it:

- **All parameters live in one concrete `NamedTuple`** (`Float64` values),
  passed to the ODE through `p` and read with literal destructuring
  `(; alpha1, …) = p`. Fully inferable, zero-cost.
- **Weather and controls are concrete structs** of callable interpolants
  (never `Dict{Symbol,Any}`).
- **Inference flexibility via a function barrier:** `update_params(base, names, x)`
  builds a new parameter `NamedTuple` with just the chosen names overridden —
  once per MCMC proposal, outside the ODE inner loop. Same freedom as globals,
  no speed penalty.

The math in `climate_rhs` and `assimilation` is a faithful, byte-for-byte port
of 2.0, so results match (see the regression check in `test/benchmark.jl`).

## Layout

    src/GreenhouseSim.jl   module + exports
    src/parameters.jl      load_params, update_params  (concrete NamedTuple)
    src/weather.jl         load_weather -> WeatherInputs
    src/controls.jl        load_measured_controls -> ControlInputs
    src/climate.jl         climate_rhs (4 climate ODEs) + psychrometrics
    src/photosynthesis.jl  assimilation (FvCB)
    src/forward_map.jl     SimBase, SimContext, rhs!, simulate, forward_map
    test/benchmark.jl      type-stability + allocation + speed + regression
    test/run_inference.jl  t-walk MCMC (mirrors 2.0 for_inference_V4)

## Run

    cd invernaderosV2.2/GreenhouseSim
    julia --project=. test/benchmark.jl        # verify speed & correctness
    julia --project=. test/run_inference.jl    # run the inference

The scripts read the config JSONs and data from the sibling `invernaderos2.0`
folder, so keep V2.2 beside it. The environment (`Project.toml` / `Manifest.toml`)
is copied from the proven 2.0 setup.

## Choosing which parameters to infer

Edit `QoI` in `test/run_inference.jl` — any names present in the JSON configs,
with a prior range in `QoI_dict`. No other change needed; `update_params`
handles the rest. Example:

    QoI = ["alpha1", "alpha4", "tau1", "eta1"]

## Next step — the growth model

`assimilation` currently returns `A` that only feeds the CO2 balance, and LAI is
a constant input (see the fidelity note in `src/photosynthesis.jl`). The daily
loop in `_run` (forward_map.jl) is the hook: integrate `A` over each day,
partition dry matter, update LAI/biomass, and feed LAI back into the next day.
