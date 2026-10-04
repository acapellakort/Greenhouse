# Invernaderos — Greenhouse Simulation & Control Monorepo

> Research project · Instituto de Matemáticas, UNAM  
> Contact: capella@im.unam.mx

## Overview

End-to-end pipeline for **physics-based greenhouse simulation**, **Bayesian crop-parameter calibration**, and **reinforcement-learning climate control**, applied to the AGC-2018 cucumber dataset.

---

## Repository layout

```
invernaderos/
├── invernaderos2.0/          # Original Python/Julia prototype (historical)
├── invernaderosV2.1/         # Intermediate Julia rewrite (superseded)
├── invernaderosV2.2/         # ★ Current simulation engine (frozen V2.2)
│   └── GreenhouseSim/        #   Julia package: ODE climate model + FvCB photosynthesis
│       ├── src/              #   Core library
│       └── test/             #   Calibration & forecast scripts
│           ├── run_inference.jl       # Climate-parameter t-walk MCMC
│           ├── calibrate_growth.jl    # Per-team crop calibration
│           ├── calibrate_all_teams.jl # Batch calibration (all 6 AGC-2018 teams)
│           └── run_forecast.jl        # Probabilistic 5-week production forecast
├── invernaderosV3.0/         # ★ Control layer (includes V2.2 as dependency)
│   └── GreenhouseControl/    #   Julia package: SAC, MPC-distilled SAC, CEM
│       ├── src/
│       └── test/
│           └── campaign_actors/       # Trained SAC policy checkpoints (.jls)
├── GreenHouse_3.0/           # Earlier Julia prototype (code only; data excluded)
├── notes/                    # Numbered research notes (Markdown)
├── docs/                     # Extended documentation
└── handoff_greenhouses_project.md   # Project context & status
```

---

## Key design decisions

| Decision | Detail |
|---|---|
| **V2.2 engine is frozen** | All new physics goes in as default-off flags; the core ODE is not modified |
| **Bayesian calibration** | t-walk MCMC via [JTwalk.jl](https://github.com/AndresMorenoP/JTwalk.jl); posteriors saved as `.jld2` files |
| **Multi-seed / split-year evaluation** | TRAIN/VAL/TEST weather years are disjoint; metrics reported over ≥3 seeds |
| **`forecast_skill`** | 0 = perfect forecast, higher = noisier; used to benchmark probabilistic forecasts |
| **SAC control** | Soft Actor-Critic agents trained against V2.2 plant model; best checkpoints in `campaign_actors/` |

---

## Quick start

### Simulation engine

```julia
cd invernaderosV2.2/GreenhouseSim
julia --project=.
julia> using Pkg; Pkg.instantiate()
julia> include("test/calibrate_growth.jl")   # ~30–60 min
julia> include("test/run_forecast.jl")        # ~5–10 min
```

### Control training

```julia
cd invernaderosV3.0/GreenhouseControl
julia --project=.
julia> include("test/train_sac.jl")
```

---

## Data

Weather data (30 historical years, 1989–2018) and AGC-2018 crop measurements are **not stored in this repository** (size). Scripts expect them under `invernaderosV2.2/GreenhouseSim/test/data/` — see `docs/` for the expected format.

---

## License

To be determined.
