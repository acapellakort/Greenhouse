# =============================================================================
# heat_audit.jl -- season heat-balance audit of the greenhouse AIR temperature
# =============================================================================
# Run:  julia --project=. test/heat_audit.jl
#
# Re-derives every term of the dT2 (air) energy balance from climate.jl at each
# solved timestep and integrates it over the season, so we see WHERE the heat
# goes (kWh/m2). Includes a consistency check: the net of all terms must equal
# the air heat capacity times the season's end-to-end temperature change (~0).

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()
include(joinpath(@__DIR__, "..", "src", "GreenhouseControl.jl"))
using .GreenhouseControl
const GS = GreenhouseControl.GreenhouseSim
using Dates, Printf, OrdinaryDiffEq

INV          = normpath(joinpath(@__DIR__, "..", "..", ".."))
climate_json = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_climate.json")
crop_json    = joinpath(INV, "invernaderos2.0", "code", "configfiles", "constants_cropphoto.json")
growth_json  = joinpath(INV, "invernaderosV2.2", "GreenhouseSim", "configfiles", "cucumber_growth.json")
meteo        = joinpath(INV, "invernaderos2.0", "code", "data", "dataset_meteo_holanda.csv")

params = load_params(climate_json, crop_json)
params = update_params(params, ["psi2","J_max"], [27800.0, 1.15e-4])
gp = load_growth_params(growth_json)
gp = update_params(gp, ["node_rate","veg_sink_max","LAI_max","set_start_day","set_rate","Wf_max"],
                       [0.091, 21.0, 2.5, 16.0, 1.5, 16.0])

start_date = DateTime(1998, 7, 11, 0, 0); ndays = 100
weather = load_weather(meteo, start_date, start_date + Day(ndays + 5); sky_from_clouds = true)
sp = daynight_setpoints(Tset_day=22+273.15, Tset_night=19+273.15, CO2_set=1200.0,
                        light_start=4.0, light_end=20.0, VentpBand=4.0, ofset=1.0, ToutMax=12+273.15)
gT = PIGains(0.6, 0.02, 5e-4); gC = PIGains(0.01, 2e-4, 5e-4)

# --- air-balance terms [W/m2] at one instant (mirrors climate.jl dT2) --------
function air_terms(T1,T2,I1,I3,I2,I4,I5,I8,T5, U1,U6,U8,U12, p)
    (; alpha2,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha12, beta3,
       epsil3,epsil4,epsil5,epsil6,eps_screen, eta1,eta2,eta3,eta6,eta7,eta8,eta10,eta11,eta15,eta16,
       gamma1, nu1,nu2,nu3,nu4,nu5,nu6,nu7,nu8, n_pipes, omega1,
       phi1,phi2, rho3, sigma, tau1,tau2,tau3, k_ground) = p
    U5 = 0.0
    n1 = nu1*(1-eta10*U5); n2 = nu3*U6; n3 = nu2*(1-eta11*U5)
    f5 = n1*n2*I8*sqrt(n3)/(2*alpha6)
    f6 = I8 < 0.25 ? 0.25*nu4 : nu4*I8
    f7t = max(omega1*nu6*(T2-I5)/(T2+I5) + n3*I8^2, 0.0)
    f7 = (U8*nu5*n1)/(2.0*alpha6)*sqrt(f7t)
    f2 = eta7 >= eta8 ? eta6*f5 : eta6*(U1*f5)
    f4 = eta6*f7
    g2 = 1 - U1*(1-tau3)
    g3 = tau1*g2*(1-0.49*pi*gamma1*phi1)*exp(-beta3*I1)
    SkyFIR = epsil3*sigma*(T2^4 - I4^4)
    r11 = epsil4*g3*epsil3*sigma*(T5^4 - I4^4)
    r12 = eps_screen*(U1*tau2)*SkyFIR
    r13 = epsil6*SkyFIR
    r10 = r11 + r12 + r13
    r8  = I2*(eta1*(tau1*eta2 + (alpha2+alpha7)*eta3) + (alpha8*eta2 + alpha9*eta3))
    h1  = 2*alpha4*I1*(T1 - T2)
    h4  = n_pipes*(1.99*pi*phi1*gamma1*abs(I3-T2)^0.32)*(I3-T2)
    h7  = rho3*alpha5*(f2 + f4 + 0.5*f6)*(T2 - I5)          # ventilation (roof f4 included; f3=0)
    h11 = 2*nu7*(T2 - T5)/(phi2 + nu8)                       # soil via nu7
    h12 = (eta15+eta16)*alpha12*U12                          # lamp heat
    Qf  = k_ground*(T2 - T5)                                 # floor via k_ground
    return (; h1, h4, h12, r8, h7, h11, r10, Qf, r11, r12, r13)
end

# --- run the season, integrate each term -------------------------------------
gs = GS.init_growth_state(gp; LAI0=0.5)
u  = copy(GreenhouseControl.U0_CONTROL)
acc = Dict(k=>0.0 for k in (:h1,:h4,:h12,:r8,:h7,:h11,:r10,:Qf,:r11,:r12,:r13))  # J/m2
airT=Float64[]; outT=Float64[]
for d in 0:ndays-1
    u[7]=0.0; wt0 = datetime2unix(start_date)+86400.0*d
    ctx = ControlContext(params, weather, sp, gs.LAI, wt0, gT, gC, params.HEAT_PIPE)
    prob = ODEProblem(GreenhouseControl.rhs_control!, u, (0.0,86400.0), ctx)
    sol = solve(prob, Tsit5(); saveat=300.0, dtmax=300.0)
    ts = sol.t; nt=length(ts)
    terms = Vector{NamedTuple}(undef, nt)
    for j in 1:nt
        t=ts[j]; td=mod(t,86400.0); tcal=t+wt0
        T1=sol[1,j];T2=sol[2,j];C1=sol[3,j];eT=sol[5,j];T5=sol[8,j]
        I2=weather.Iglobal(tcal);I4=weather.Tsky(tcal);I5=weather.Tout(tcal);I8=weather.WindSpeed(tcal)
        Tset=sp.Tset(td);pBand=sp.VentpBand(td);ofs=sp.ofset(td);ToutMax=sp.ToutMax(td)
        ScreenRad=sp.ScreenRad(td);VPDmin=sp.VPDmin(td);VPD=(GS.Pws(T2)-sol[4,j])/1000.0
        U1 = screen_control(I2,I5,VPD,ScreenRad,ToutMax,VPDmin)
        U6,U8 = vent_control(Tset,T2,pBand,ofs)
        Uhv = humidity_vent(VPD,VPDmin); U6=max(U6,Uhv); U8=max(U8,Uhv)
        U12 = sp.Light_on(td)
        Uheat = clamp(gT.Kp*(Tset-T2)+gT.Ki*eT,0.0,1.0); I3=Uheat*params.HEAT_PIPE+(1-Uheat)*T2
        terms[j] = air_terms(T1,T2,gs.LAI,I3,I2,I4,I5,I8,T5,U1,U6,U8,U12,params)
        push!(airT,T2-273.15); push!(outT,I5-273.15)
    end
    for j in 1:nt-1
        dt=ts[j+1]-ts[j]
        for k in keys(acc); acc[k]+=0.5*(getproperty(terms[j],k)+getproperty(terms[j+1],k))*dt; end
    end
    # advance crop with this day's Pg (reuse assimilation integral)
    Aint=sol[7,end]-sol[7,1]; Pgd=Aint*1e-3*(30/44); Tm=sum(@view sol[2,:])/nt-273.15
    GS.grow!(gs,Pgd,Tm,gp,d); global u=sol[:,end]
end

kWh(x)= x/3.6e6
println("\n=== season air heat-balance audit (1998, 100 d) ===")
println("mean air ", round(sum(airT)/length(airT),digits=1), " C   mean outside ", round(sum(outT)/length(outT),digits=1), " C")
println("\nGAINS  [kWh/m2]")
@printf("  pipe heating h4   %8.1f\n", kWh(acc[:h4]))
@printf("  solar to air r8   %8.1f\n", kWh(acc[:r8]))
@printf("  lamp heat  h12    %8.1f\n", kWh(acc[:h12]))
@printf("  canopy->air h1    %8.1f\n", kWh(acc[:h1]))
println("LOSSES [kWh/m2]")
@printf("  ventilation h7    %8.1f\n", kWh(acc[:h7]))
@printf("  FIR to sky  r10   %8.1f   (r11 floor %.1f | r12 screen %.1f | r13 cover-air %.1f)\n", kWh(acc[:r10]), kWh(acc[:r11]), kWh(acc[:r12]), kWh(acc[:r13]))
@printf("  soil (nu7)  h11   %8.1f\n", kWh(acc[:h11]))
@printf("  floor (k_ground)  %8.1f\n", kWh(acc[:Qf]))
net = acc[:h1]+acc[:h4]+acc[:h12]+acc[:r8]-acc[:h7]-acc[:r10]-acc[:h11]-acc[:Qf]
aircap = params.phi2*params.rho3*params.alpha5
@printf("\nNET (gains-losses) = %.1f kWh/m2   [consistency: air_cap*dT2/season ~ %.3f kWh/m2]\n",
        kWh(net), kWh(aircap*(airT[end]-airT[1])))
