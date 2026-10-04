# =============================================================================
# run_env_test.jl -- validate the stepping RL env against season_economics
# =============================================================================
# Run:  julia --project=. test/run_env_test.jl
# Replays the fixed baseline policy through env_step! and checks the summed daily
# rewards equal the whole-season net profit (env is consistent).

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
using Dates

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params = update_params(load_params(climate_json, crop_json), ["psi2","J_max"], [27800.0, 1.15e-4])
gp = update_params(load_growth_params(growth_json),
        ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max"],
        [0.091, 21.0, 2.5, 16.0, 1.5, 16.0])

start_date = DateTime(1998, 7, 11, 0, 0); ndays = 100
weather = load_weather(meteo, start_date, start_date + Day(ndays + 5); sky_from_clouds = true)

env = TwinEnv(params, gp, weather, datetime2unix(start_date), ndays)

# baseline policy as a normalized action (Tday22,Tnight19,CO2 1200,lh16,pBand4,ofs1,ToutMax12)
base_phys = [22.0, 19.0, 1200.0, 16.0, 4.0, 1.0, 12.0, 0.5, 50.0, 4.0, 1.0]  # +VPDmin,ScreenRad,light_start,lamp_intensity
a = (base_phys .- ACT_LO) ./ (ACT_HI .- ACT_LO)
println("baseline action (normalized): ", round.(a, digits=3))

obs = env_reset!(env)
total = 0.0; done = false
Rheat=0.0; Rco2=0.0; Relec=0.0
while !done
    obs, r, done, info = env_step!(env, a)
    global total += r
    global Rheat += info.hkWh; global Rco2 += info.ckg; global Relec += info.elec
end
println("\n--- env rollout (baseline policy) ---")
println("summed daily reward = ", round(total, digits=2), " EUR/m2   (should match run_baseline NET ~ -33.2)")
println("final yield = ", round(env.gs.yield_FW, digits=2), " kg/m2")
println("heating ", round(Rheat,digits=1), " kWh | CO2 ", round(Rco2,digits=1), " kg | elec ", round(Relec,digits=1), " kWh")
println("obs (final) = ", round.(obs, digits=3))
