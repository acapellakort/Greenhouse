# Performance Assessment — V2.1 forward map (speed + inference flexibility)

_Question: is V2.1 heading the right way to make the forward-map simulations faster
and flexible for parameter inference, while keeping the ability to infer any
parameter by name (2.0 did this with global variables)?_

**Short answer:** the *packaging* and the *`merge`-based parameter update* are the
right direction and worth keeping. But V2.1 does **not** yet make the forward map
faster, because the physics kernel still reads non-`const` global variables — the
same performance killer as 2.0. The "infer any parameter" goal is fully preservable
**without** globals, and globals are precisely what blocks the speedup.

---

## 1. Decisive finding

`invernaderosV2.1/GreenhouseSim/src/utils/import_constants_NoUnits.jl` still defines
**all** physics constants as non-`const` globals via:

```julia
eval(:(global $global_var_name = $unitless_value))
```

`rhs_fast` (called ~10^5–10^6 times per solve) reads ~60 of these as free variables:
`gamma5, delta1, delta4, delta5, rho3, alpha5, nu1, eta10, phi7, gamma2, gamma3, …`.
Only `eta1, tau1, eta2` come from the new `ModelParams` struct.

⇒ The inner loop is essentially identical to 2.0, so **runtime is essentially
unchanged.** This is almost certainly the source of the README's "13.33 GiB
allocations" for the full model — a signature of type instability from globals.

## 2. Why non-`const` globals are slow in Julia

A non-`const` global's type may change at any time, so the compiler must assume
`Any`: every read is boxed, every operation is dynamically dispatched, and no
inlining/SIMD/optimization survives. In a function called hundreds of thousands of
times per simulation this dominates runtime and generates huge heap allocation.
Removing this is the single biggest speed lever available.

## 3. What V2.1 got RIGHT (keep these)

- Proper package layout (`Project.toml`, module, `SimulationContext`).
- Removed the fatal `eval(:(global …))` **inside the MCMC loop** (2.0's `FM_SimLoop!`
  set globals per proposal via `eval` in global scope → recompilation + world-age
  cost every step). V2.1 replaces it with `merge` of a parameter object — correct.
- `FM_SimLoop!` already sketches the right idea:
  `merge(base_ctx.params.vals, NamedTuple{Tuple(Symbol.(theta_names))}(x))`.

## 4. What V2.1 got WRONG / left type-unstable

- **Physics still reads globals** (see §1) — the wrapper struct is bypassed by
  `rhs_fast`. The struct helps speed only if the RHS reads *all* params from it.
- `struct ModelParams; vals::NamedTuple; …` — bare `NamedTuple` is an **abstract**
  type. Field access is not inferrable. Needs a type parameter, e.g.
  `struct ModelParams{NT<:NamedTuple}; vals::NT; …end` — or just use a bare concrete
  NamedTuple and drop the wrapper.
- `Base.getproperty(::ModelParams, s::Symbol) = p.vals[s]` — indexing a NamedTuple by
  a **runtime** Symbol is type-unstable. Fine outside hot loops; not inside the RHS.
- `SimulationContext.weather::Dict{Symbol, Any}` — `Any` value type ⇒ every
  `weather[:Tout](t)` in the RHS returns `Any` ⇒ unstable on the hot path. Use a
  concrete struct or a typed NamedTuple of interpolant objects.
- V2.1 also does not currently run (missing `load_weather_data`, missing data files,
  `FM_SimLoop!` has no return) — see note 01.

## 5. The design that is BOTH fast AND flexible

**Requirement to preserve:** choose *any* subset of parameters to infer, by name, at runtime.

1. **One concrete parameter container**, passed through the ODE `p`:
   a flat `NamedTuple` of ~70 `Float64` values (fully concrete, composes with `merge`).
2. **RHS reads via literal `@unpack`:** `@unpack alpha1, alpha2, … = p` → type-stable, zero-cost.
   No globals anywhere in `rhs_fast` / `climate_model_rhs` / photosynthesis.
3. **Inference flexibility via `merge` (function barrier):**
   ```julia
   p = merge(base_params, NamedTuple{Tuple(Symbol.(QoI))}(Tuple(x)))
   solve(remake(prob; p=p, u0=u0), Tsit5(), saveat=t_eval)
   ```
   Overrides exactly the inferred params, chosen at runtime by name — the same
   freedom globals gave, but the (tiny) dynamic cost is once per MCMC step, outside
   the ODE inner loop. This keeps your key feature with no speed penalty.
4. **Weather & controls as concrete typed containers** (struct or typed NamedTuple of
   interpolation objects), not `Dict{…,Any}`.

## 6. Secondary optimizations (after §5)

- `minimum([Ar_j, Ar_c, Ar_p])` → `min(Ar_j, Ar_c, Ar_p)` (removes a heap alloc per
  photosynthesis call); same for any `min([...])`/`max([...])`.
- Don't build `controls = [...]` / `inputs = [...]` arrays each RHS call (2.0 & 3.0 do
  this) — pass tuples/scalars. V2.1's `physics_ode.jl` already passes scalars ✔.
- Daily loop: `remake(prob; u0=…, p=…)` instead of constructing a new `ODEProblem`
  each day; reuse the compiled problem across MCMC steps.
- Ensure the RHS is fully in-place (`du` ✔) and closures capture no non-`const` globals.
- Re-benchmark solver choice once type-stable (Tsit5 kept for convergence today).
- For inference throughput: parallelize chains / forward solves once per-solve is fast.

## 7. How to verify

- `@code_warntype rhs_fast(…)` → no red `Any`/`Union` on hot variables.
- `@btime` a single RHS call (**target: 0 allocations**) and a full-day `solve`.
- `julia --track-allocation=user` on a full sim → allocations should drop from GB to KB–MB.
- Regression check: forward map on fixed params must match 2.0 output within solver tol.

## 8. Recommendation

Keep V2.1's structure and its `merge` inference pattern. Finish it by routing **all**
parameters (and weather/controls) through concrete typed containers into the RHS, and
**delete `import_constants_NoUnits.jl`'s global importer.** Expect ~10–100× on the
forward map — the leverage inference and agent training need. Open decision: refactor
V2.1 in place, or branch a clean V2.2 built from 2.0's proven physics.
