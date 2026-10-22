# =============================================================================
# climate_calibration.jl -- likelihood of the n79 climate calibration (port of run_climate_calibration_cl.jl)
# =============================================================================
# Design: 7-day blocks with one unscored spin-up day and the measured initial state.
#   TRAIN 7-13, 28-34, 63-69, 91-97, 110-114 | VAL 42-48, 77-83 | TEST the rest | COMMON TEST 49-62, 70-76, 84-90, 98-109
# Observations: hourly-bin means of air T and RH (bins with >= half of their points observed), and dark-hour CO2
# (whole hour dark, lamps off, valve shut). Energy = sum_k [n_k log sigma_k + SSE_k / (2 sigma_k^2)].

"Free climate parameters of n79 (order of the V2.2 calibration), followed by the noise terms."
const CLIMATE_QOI = ["gamma4", "vent_scr_perm", "k_ground", "nu1", "nu2", "C_soil", "k_deep", "C_mass", "h_mass", "f_mass",
                     "lamp_fcan", "vs_gwind", "vs_cwind", "vs_exp", "nu4", "h_cover_out", "air_cap_mult", "leak_cs"]
const CLIMATE_NOISE = ["sigma_T", "sigma_RH", "sigma_C"]
const TRAIN_BLOCKS = [(7, 13), (28, 34), (63, 69), (91, 97), (110, 114)]
const VAL_BLOCKS   = [(42, 48), (77, 83)]
const NPD = 289                                         # 5-min points per day incl. both ends
const TEVAL = 0.0:300.0:86400.0

struct Block
    days::Vector{Float64}; u0::Vector{Float64}
    binsT::Vector{Vector{Int}}; obsT::Vector{Float64}
    binsR::Vector{Vector{Int}}; obsR::Vector{Float64}
    binsC::Vector{Vector{Int}}; obsC::Vector{Float64}
end

_binmean(v, ix) = (s = 0.0; @inbounds for i in ix; s += v[i]; end; s / length(ix))

"Block (a, b): days a-1 (spin-up) .. b, measured initial state at the start of the spin-up day."
function make_block(S, a, b)
    A(v) = StepFun(S.ts, v)
    AT = A(S.Tair_obs); AR = A(S.RH_obs); AC = A(S.CO2_obs)
    t0 = S.t0
    d0 = max(a - 1, 0); ds = collect(Float64, d0:b)
    tstart = t0 + 86400.0 * d0
    Ta = AT(tstart); Ta = isnan(Ta) ? S.Tair_obs[findfirst(!isnan, S.Tair_obs)] : Ta
    Rh = AR(tstart); Rh = isnan(Rh) ? 80.0 : Rh
    Co = AC(tstart); Co = isnan(Co) ? 575.0 : Co
    ts_ = filter(!isnan, AT.(tstart .+ (0.0:300.0:86400.0))); Td0 = isempty(ts_) ? Ta : mean(ts_)
    u0b = [Ta, Ta, Co, Rh / 100 * Pws(Ta), 0.0, Ta, Td0, Td0]
    times = vcat(tstart, [t0 + d * 86400.0 + tj for d in ds for tj in TEVAL])
    To = AT.(times); Ro = AR.(times); Co_ = AC.(times)
    darkC = (A(S.Ig).(times) .< 1.0) .& (A(S.U12).(times) .< 0.01) .& (A(S.U10).(times) .<= 0.0)
    bT = Vector{Int}[]; oT = Float64[]; bR = Vector{Int}[]; oR = Float64[]; bC = Vector{Int}[]; oC = Float64[]
    for (k, d) in enumerate(ds)
        Int(d) < a && continue                                   # spin-up day not scored
        for hb in 0:23
            cols = [2 + (k - 1) * NPD + j for j in 0:NPD-1 if min(floor(Int, j * 300 / 3600), 23) == hb]
            okT = [c for c in cols if !isnan(To[c])]; okR = [c for c in cols if !isnan(Ro[c])]
            length(okT) >= length(cols) / 2 && (push!(bT, okT); push!(oT, mean(To[okT])))
            length(okR) >= length(cols) / 2 && (push!(bR, okR); push!(oR, mean(Ro[okR])))
            okC = [c for c in cols if !isnan(Co_[c]) && darkC[c]]
            length(okC) == length(cols) && (push!(bC, okC); push!(oC, mean(Co_[okC])))
        end
    end
    return Block(ds, u0b, bT, oT, bR, oR, bC, oC)
end

"TRAIN blocks of the n79 design."
climate_blocks(S; blocks = TRAIN_BLOCKS) = [make_block(S, a, b) for (a, b) in blocks]

"(SSE_T, SSE_RH, SSE_CO2) of one block, Inf if the solve fails."
function block_sse(p, S, b::Block, rtol)
    us = simdays(p, S, b.days, b.u0, rtol); us === nothing && return (Inf, Inf, Inf)
    T2 = us[2, :]; RH = rhf.(us[2, :], us[4, :]); C1 = us[3, :]
    sT = 0.0; sR = 0.0; sC = 0.0
    for k in eachindex(b.binsC); sC += (_binmean(C1, b.binsC[k]) - b.obsC[k])^2; end
    for k in eachindex(b.binsT); sT += (_binmean(T2, b.binsT[k]) - b.obsT[k])^2; end
    for k in eachindex(b.binsR); sR += (_binmean(RH, b.binsR[k]) - b.obsR[k])^2; end
    f(s) = isfinite(s) ? s : Inf
    return (f(sT), f(sR), f(sC))
end

"""
    climate_energy(z, p0, S, blocks; rtol) -> U

z = values of CLIMATE_QOI followed by sigma_T, sigma_RH, sigma_C; p0 = base parameters (load_params of
climate_n79.json); S = load_agc_season(...); blocks = climate_blocks(S).
"""
function climate_energy(z, p0, S, blocks; rtol = 1e-6)
    nm = length(CLIMATE_QOI)
    p = update_params(p0, CLIMATE_QOI, z[1:nm])
    sT = 0.0; sR = 0.0; sC = 0.0
    for b in blocks
        a, c, g = block_sse(p, S, b, rtol); (isfinite(a) && isfinite(c) && isfinite(g)) || return Inf
        sT += a; sR += c; sC += g
    end
    nT = sum(length(b.obsT) for b in blocks); nR = sum(length(b.obsR) for b in blocks); nC = sum(length(b.obsC) for b in blocks)
    σT, σR, σC = z[nm+1], z[nm+2], z[nm+3]
    return nT * log(σT) + sT / (2σT^2) + nR * log(σR) + sR / (2σR^2) + nC * log(σC) + sC / (2σC^2)
end

# --- full-season evaluation ------------------------------------------------------------------------------------------
const COMMON_TEST = vcat(collect(49:62), collect(70:76), collect(84:90), collect(98:109))

"Day sets of the n79 design: TRAIN, VAL, TEST, COMMON TEST."
function design_sets(ndays = 115)
    tr = reduce(vcat, [collect(a:b) for (a, b) in TRAIN_BLOCKS]); va = reduce(vcat, [collect(a:b) for (a, b) in VAL_BLOCKS])
    return [("TRAIN", tr), ("VAL", va), ("TEST", setdiff(0:ndays-1, tr, va)), ("COMMON TEST", COMMON_TEST)]
end

"""
    season_eval(p, S; rtol) -> (us, table::DataFrame)

Simulate the whole season from S.u0 and compare with the observations on the 5-min grid: air T (RMSE, bias, night,
day), vapour pressure bias (day, night, hPa), RH bias and CO2 RMSE / bias, per day set. Day: Iglob > 50 W m-2; night: < 1.
"""
function season_eval(p, S; rtol = 1e-6)
    days = collect(0.0:(S.ndays - 1))
    us = simdays(p, S, days, S.u0, rtol); us === nothing && error("full-season solve failed")
    A(v) = StepFun(S.ts, v)
    pt = vcat(S.t0, [S.t0 + d * 86400.0 + tj for d in days for tj in TEVAL])
    To = A(S.Tair_obs).(pt); Ro = A(S.RH_obs).(pt); Co = A(S.CO2_obs).(pt); Ig = A(S.Ig).(pt)
    cd = vcat(-1, [(c - 2) ÷ NPD for c in 2:length(pt)])
    T2 = us[2, :]; RH = rhf.(us[2, :], us[4, :]); C1 = us[3, :]
    VPm = RH ./ 100 .* Pws.(T2) ./ 100; VPo = Ro ./ 100 .* Pws.(To) ./ 100
    nite = Ig .< 1.0; dayt = Ig .> 50.0
    rows = []
    mb(v, o, k) = any(k) ? mean(v[k] .- o[k]) : NaN
    for (nm, ds) in design_sets(S.ndays)
        ins = in.(cd, Ref(Set(ds))); okT = .!isnan.(To) .& ins; okT[1] = false
        okR = okT .& .!isnan.(Ro); okC = .!isnan.(Co) .& ins; okC[1] = false
        push!(rows, (set = nm, T_RMSE = sqrt(mean((T2[okT] .- To[okT]) .^ 2)), T_bias = mb(T2, To, okT),
                     T_night = mb(T2, To, okT .& nite), T_day = mb(T2, To, okT .& dayt),
                     VP_day_hPa = mb(VPm, VPo, okR .& dayt), VP_night_hPa = mb(VPm, VPo, okR .& nite),
                     RH_bias = mb(RH, Ro, okR), CO2_RMSE = sqrt(mean((C1[okC] .- Co[okC]) .^ 2)), CO2_bias = mb(C1, Co, okC)))
    end
    return us, DataFrame(rows)
end
