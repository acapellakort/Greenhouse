# =============================================================================
# rl_env.jl  --  stepping RL environment around the digital twin (Stage 3b)
# =============================================================================
# Gym-style: env_reset!(env) -> obs; env_step!(env, action) -> (obs, reward, done, info).
# One step = one day. Action = 9 daily setpoints (normalized [0,1]); density is an
# episode-level schedule set at reset. Reward = daily net profit EUR/m2. Each episode
# samples a covered weather year (WeatherBank) and the obs includes a noisy forecast.
#
# Crop model is switchable:
#   :bigleaf -> GrowthState/grow!, daily Pg from the climate ODE (state 7).
#   :cohort  -> CohortState/grow_cohort!, daily Pg from the LAYERED per-cohort
#               canopy on the day's simulated climate (note 12; ~realistic yield).

using OrdinaryDiffEq, Random

const ACT_LO    = [16.0, 14.0, 400.0,  0.0,  2.0, 0.5,  5.0, 0.2,   5.0,  0.0, 0.0]
const ACT_HI    = [30.0, 22.0, 1200.0, 20.0, 12.0, 4.0, 22.0, 1.2, 100.0, 12.0, 1.0]
const ACT_NAMES = ["Tset_day","Tset_night","CO2_set","light_hours","VentpBand",
                   "ofset","ToutMax","VPDmin","ScreenRad","light_start","lamp_intensity"]

function decode_setpoints(a)
    v = ACT_LO .+ clamp.(a, 0.0, 1.0) .* (ACT_HI .- ACT_LO)
    Tday, Tnight, CO2, lh, pBand, ofs, ToutMaxC, VPDmin, ScreenRad, ls, lint = v
    return daynight_setpoints(Tset_day = Tday + 273.15, Tset_night = Tnight + 273.15,
                              CO2_set = CO2, light_start = ls, light_end = clamp(ls + lh, ls, 24.0),
                              VentpBand = pBand, ofset = ofs, ToutMax = ToutMaxC + 273.15,
                              ScreenRad = ScreenRad, VPDmin = VPDmin, lamp_intensity = lint)
end

# --- environment ----------------------------------------------------------------
mutable struct TwinEnv
    params::Any; gp::Any; weather::Any
    start_unix::Float64; ndays::Int
    gT::PIGains; gC::PIGains; LAI0::Float64
    crop::Symbol
    gs::Any; u::Vector{Float64}; day::Int; prev_yield::Float64; prev_dens::Float64
    density_sched::Function
    bank::Any; forecast::Any; rng::Any; forecast_skill::Float64; year::Int
    last_T::Float64; last_CO2::Float64; last_VPD::Float64
end

const _FC_H = 3

_init_crop(gp, crop, LAI0) = crop === :cohort ? init_cohort_state(gp; LAI0=LAI0) :
                                                init_growth_state(gp; LAI0=LAI0)
_wleaf(gs, gp) = hasproperty(gs, :W_leaf) ? gs.W_leaf : gs.LAI / gp.SLA

function _newenv(params, gp, weather, start_unix, ndays, bank, crop; gT, gC, LAI0, forecast_skill, rng)
    TwinEnv(params, gp, weather, start_unix, ndays, gT, gC, LAI0, crop,
            _init_crop(gp, crop, LAI0), copy(U0_CONTROL), 0, 0.0, gp.stem_density,
            (d -> gp.stem_density), bank, nothing, rng, forecast_skill, 0,
            20.0, 800.0, 0.8)
end

TwinEnv(params, gp, weather, start_unix::Float64, ndays::Int;
        gT = PIGains(0.6, 0.02, 5e-4), gC = PIGains(0.01, 2e-4, 5e-4), LAI0 = 0.5,
        crop = :bigleaf, forecast_skill = 1.0, rng = Random.default_rng()) =
    _newenv(params, gp, weather, start_unix, ndays, nothing, crop;
            gT = gT, gC = gC, LAI0 = LAI0, forecast_skill = forecast_skill, rng = rng)

TwinEnv(params, gp, bank::WeatherBank;
        ndays = bank.ndays, gT = PIGains(0.6, 0.02, 5e-4), gC = PIGains(0.01, 2e-4, 5e-4),
        LAI0 = 0.5, crop = :bigleaf, forecast_skill = 1.0, rng = Random.default_rng()) =
    _newenv(params, gp, nothing, 0.0, ndays, bank, crop;
            gT = gT, gC = gC, LAI0 = LAI0, forecast_skill = forecast_skill, rng = rng)

function observe(env::TwinEnv)
    gs = env.gs
    crop = Float64[ env.day / env.ndays, gs.LAI, gs.yield_FW, gs.node, _wleaf(gs, env.gp),
                    isempty(gs.fruits) ? 0.0 : sum(f.W for f in gs.fruits) ]
    clim = Float64[ env.last_T / 30.0, env.last_CO2 / 1200.0, env.last_VPD / 1.5 ]
    fc = env.forecast === nothing ? zeros(_FC_H * 3) : forecast_row(env.forecast, env.day)
    fcn = similar(fc)
    @inbounds for l in 0:(_FC_H - 1)
        fcn[3l+1] = fc[3l+1] / 30.0; fcn[3l+2] = fc[3l+2] / 400.0; fcn[3l+3] = fc[3l+3] / 10.0
    end
    return vcat(crop, clim, fcn)
end

function env_reset!(env::TwinEnv; density_sched = nothing, year = nothing)
    if env.bank !== nothing
        w, su, y = sample_window(env.bank; rng = env.rng, year = year)
        env.weather = w; env.start_unix = su; env.year = y
    end
    env.forecast = SeasonForecast(env.weather, env.start_unix, env.ndays;
                                  H = _FC_H, skill = env.forecast_skill, rng = env.rng)
    env.density_sched = density_sched === nothing ? (d -> env.gp.stem_density) : density_sched
    dens0 = env.density_sched(0)
    env.gs = _init_crop(env.gp, env.crop, env.LAI0 * dens0 / 2.5)
    env.u  = copy(U0_CONTROL); env.day = 0; env.prev_yield = 0.0; env.prev_dens = dens0
    env.last_T = 20.0; env.last_CO2 = 800.0; env.last_VPD = 0.8
    return observe(env)
end

"Hourly canopy climate matrix [Tair_K; CO2_mg/m3; Iw; VPD_kPa] from the day's ODE solution."
function _pg_matrix(sol, weather, sp, wt0, p)
    (; eta1, tau1, eta2, alpha12) = p
    trans = (1 - eta1) * tau1 * eta2
    nt = length(sol.t); M = Matrix{Float64}(undef, 4, nt)
    @inbounds for j in 1:nt
        tabs = wt0 + sol.t[j]; td = mod(sol.t[j], 86400.0)
        M[1, j] = sol[2, j]                                   # air temperature [K]
        M[2, j] = sol[3, j]                                   # CO2 [mg/m3]
        M[3, j] = trans * weather.Iglobal(tabs) + alpha12 * sp.Light_on(td)   # canopy light
        M[4, j] = (Pws(sol[2, j]) - sol[4, j]) / 1000.0        # leaf-air VPD [kPa]
    end
    return M
end

function env_step!(env::TwinEnv, action)
    sp = decode_setpoints(action)
    dens = env.density_sched(env.day)
    # management thinning (cohort only): a density drop pulls stems (leaves + fruit)
    if env.crop === :cohort && dens < env.prev_dens - 1e-9
        r = dens / env.prev_dens
        for c in env.gs.leaves; c.area *= r; end
        for f in env.gs.fruits; f.n    *= r; end
        env.gs.LAI *= r
    end
    env.u[7] = 0.0
    wt0 = env.start_unix + 86400.0 * env.day
    ctx = ControlContext(env.params, env.weather, sp, env.gs.LAI, wt0, env.gT, env.gC, env.params.HEAT_PIPE)
    prob = ODEProblem(rhs_control!, env.u, (0.0, 86400.0), ctx)
    sol  = solve(prob, Tsit5(); saveat = 3600.0, dtmax = 300.0)
    _bad(so) = so.t[end] < 86399.0 || !all(isfinite, @view so[:, end])
    if _bad(sol)                                   # stiff-solver fallback before penalizing
        sol = solve(prob, Rosenbrock23(); saveat = 3600.0, dtmax = 300.0)
    end
    # robustness guard: an unstable solve (blow-up / dt->eps abort) leaves an
    # incomplete or non-finite trajectory. Penalize, keep the last good state,
    # skip growth, and continue -- so training/eval never ingest garbage.
    if _bad(sol)
        env.day += 1
        done = env.day >= env.ndays
        return observe(env), -3.0, done,
               (; hkWh = 0.0, ckg = 0.0, elec = 0.0, dyield = 0.0, LAI = env.gs.LAI,
                  Tmean = env.last_T, CO2mean = env.last_CO2, VPDmean = env.last_VPD,
                  density = dens, year = env.year, Pg = 0.0, failed = true)
    end
    n = size(sol, 2)
    Tm   = sum(@view sol[2, :]) / n - 273.15
    Cm   = sum(@view sol[3, :]) / n
    VPDm = sum((Pws(sol[2, j]) - sol[4, j]) / 1000.0 for j in 1:n) / n
    hkWh, ckg, lpk, loff = day_resources(sol, sp, env.params, env.gT, env.gC)

    gpd = update_params(env.gp, ("stem_density",), (dens,))
    if env.crop === :cohort
        M   = _pg_matrix(sol, env.weather, sp, wt0, env.params)
        Pgd = daily_assimilation_layered(M, env.gs.leaves, env.params, 3600.0, env.gp.k_ext)
        q10e = env.gp.resp_diurnal == 0.0 ? NaN :
               sum(env.gp.Q10_resp ^ ((M[1, j] - 273.15 - 25.0) / 10.0) for j in 1:size(M, 2)) / size(M, 2)
        grow_cohort!(env.gs, Pgd, Tm, gpd, env.day; dens = dens / 2.5, q10_resp = q10e)
    else
        Pgd = (sol[7, end] - sol[7, 1]) * 1e-3 * (30.0 / 44.0)
        grow!(env.gs, Pgd, Tm, gpd, env.day)
    end
    env.u = sol[:, end]; env.prev_dens = dens
    env.last_T = Tm; env.last_CO2 = Cm; env.last_VPD = VPDm

    dyield = env.gs.yield_FW - env.prev_yield; env.prev_yield = env.gs.yield_FW
    reward = max(dyield, 0.0) * PRICE_CUKE - hkWh * COST_HEAT - ckg * COST_CO2 -
             (lpk * COST_ELEC_PEAK + loff * COST_ELEC_OFF) - FIXED_GH_PER_YR / 365.0 -
             (env.day == 0 ? PLANT_COST_STEM * dens : 0.0)

    env.day += 1
    done = env.day >= env.ndays
    return observe(env), reward, done,
           (; hkWh, ckg, elec = lpk + loff, dyield, LAI = env.gs.LAI, Tmean = Tm,
              CO2mean = Cm, VPDmean = VPDm, density = dens, year = env.year, Pg = Pgd)
end
