# =============================================================================
# run_control.jl  --  Stage-1 demo: twin tracking a programmed climate
# =============================================================================
# Run:  julia --project=. test/run_control.jl
#
# Simulates a few days under a day/night setpoint program and plots how the
# controlled air temperature and CO2 track their setpoints, PLUS the actuator
# signals (heating-pipe temperature, vent opening, CO2 dosing valve) so we can
# see WHY tracking succeeds or fails and tune the controllers accordingly.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Plots, Dates

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params     = load_params(climate_json, crop_json)
# Equipment/physics reconciliation for the CONTROL twin (override, keeps V2.2 frozen):
#  - psi2: CO2 dosing capacity = 100 kg/(ha*h) = 2.78 mg/m2/s (was ~2x undersized).
#  - nu4 / nu4CO2: heat/CO2 LEAKAGE were inflated ~100x (0.01 vs Vanthoor's ~1e-4) to
#    stabilize the ODE in the fast forward map. That artificial leak is the shared cause
#    of the night-temperature undershoot AND the CO2 that won't enrich. Reduce 10x here.
params     = update_params(params,
                ["psi2", "nu4", "nu4CO2"],
                [27800.0, 1.0e-4, 1.0e-4])
start_date = DateTime(1998, 7, 11, 0, 0)
ndays      = 5
weather    = load_weather(meteo, start_date, start_date + Day(ndays + 1); sky_from_clouds = true)

# --- setpoint program (agent action surface). CO2 setpoint is ABOVE ambient
#     (outside ~835 mg/m3) so enrichment/dosing actually engages. -------------
sp = daynight_setpoints(
    Tset_day = 21 + 273.15, Tset_night = 18 + 273.15,
    CO2_set  = 1500.0,
    light_start = 6.0, light_end = 20.0,
    VentpBand = 4.0, ofset = 0.5, ToutMax = 22 + 273.15,
)

# --- controller gains (tunable) ----------------------------------------------
gT = PIGains(0.6, 0.02, 0.0005)     # heating: stronger than the first guess
gC = PIGains(0.01, 2e-4, 0.0005)    # CO2 dosing
HEAT_PIPE = params.HEAT_PIPE        # 363.15 K (90 C)

sol = run_control(params, weather, sp, datetime2unix(start_date), ndays; lai = 2.0, gT = gT, gC = gC)

# --- recompute actuators post-hoc for diagnostics ----------------------------
t     = sol.t
hours = t ./ 3600
T2    = [u[2] for u in sol.u] .- 273.15
C1    = [u[3] for u in sol.u]
Tset  = [sp.Tset(mod(tt, 86400.0)) for tt in t] .- 273.15
CO2s  = [sp.CO2_set(mod(tt, 86400.0)) for tt in t]

Ipipe = Float64[]; Vent = Float64[]; Dose = Float64[]
for (tt, u) in zip(t, sol.u)
    td = mod(tt, 86400.0)
    T2K = u[2]; C1v = u[3]; eT = u[5]; eC = u[6]
    Tk = sp.Tset(td)
    errT = Tk - T2K
    Uh = clamp(gT.Kp*errT + gT.Ki*eT, 0.0, 1.0)
    push!(Ipipe, Uh*HEAT_PIPE + (1-Uh)*T2K - 273.15)               # pipe temp [C]
    us, ur = vent_control(Tk, T2K, sp.VentpBand(td), sp.ofset(td))
    push!(Vent, 100*(us + ur)/2)                                   # mean vent %
    errC = sp.CO2_set(td) - C1v
    push!(Dose, 100*clamp(gC.Kp*errC + gC.Ki*eC, 0.0, 1.0))        # dosing %
end

figdir = joinpath(@__DIR__, "figures"); mkpath(figdir)
p1 = plot(hours, T2, lw=2, label="T_air", ylabel="deg C", title="Air temperature vs setpoint")
plot!(p1, hours, Tset, lw=2, ls=:dash, label="Tset")
p2 = plot(hours, C1, lw=2, label="CO2 air", ylabel="mg/m3", title="CO2 vs setpoint")
plot!(p2, hours, CO2s, lw=2, ls=:dash, label="CO2 set")
p3 = plot(hours, Ipipe, lw=2, label="pipe temp (C)", xlabel="hour", ylabel="signal",
          title="Actuators")
plot!(p3, hours, Vent, lw=2, label="vent %")
plot!(p3, hours, Dose, lw=2, label="CO2 dosing %")
savefig(plot(p1, p2, p3, layout=(3,1), size=(950, 950)), joinpath(figdir, "control_tracking.png"))

println("mean |T2 - Tset|   = ", round(sum(abs.(T2 .- Tset))/length(T2), digits=2), " C")
println("mean |C1 - CO2set| = ", round(sum(abs.(C1 .- CO2s))/length(C1), digits=1), " mg/m3")
println("mean pipe temp = ", round(sum(Ipipe)/length(Ipipe), digits=1), " C ; ",
        "max pipe = ", round(maximum(Ipipe), digits=1), " C")
println("Figure: ", joinpath(figdir, "control_tracking.png"))
