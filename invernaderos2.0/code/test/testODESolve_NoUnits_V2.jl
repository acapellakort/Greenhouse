using CSV
using DataFrames
using Dates
using Unitful
using DifferentialEquations
using Distributed
using BenchmarkTools
using LaTeXStrings




# Include the script that imports constants
include("../utils/import_constants_NoUnits.jl")

# include all functions
include("../utils/climate_functionsV2_NoUnits.jl")
include("../utils/cropphoto_functions_NoUnits.jl")

# Import weather data


include("../utils/readweatherdata_NoUnits.jl")

include("../utils/temp_controls.jl")

include("../utils/control_functions.jl")

include("../utils/import_control_set.jl")

# Define the right-hand side function of the ODE
#@everywhere 
function rhsPLL_PID!(du, u, p, t)

    # Unpack  y = [T1, T2, C1, V1] and u[5] = integral_error 
    # 84.715 ns
    T1 = u[1]
    T2 = u[2] 
    C1 = u[3] 
    V1 = u[4] 
    integral_error =  u[5]   

    # from simulation time t to date(unix) time t_calendar #18.784 ns
    t_calendar = t  + p  

    #weather data is given in calendar days   12.58 ns
     I1 =   LAIf(t_calendar)
     I2 =   Iglobalf(t_calendar)  
     I4 =   Tskyf(t_calendar) 
     I5 =   Toutf(t_calendar) 
     I6 =   TmechCoolf(t_calendar) 
     I7 =   Tsoilf(t_calendar) 
     I8 =   WindSpeedf(t_calendar) 
     I9 =   (1 - eta1) * tau1 * eta2 * Idocelf(t_calendar)  
     I10 =  CO2outf(t_calendar)
     I11 =  VPoutf(t_calendar)
    
     #Assimilates function photosynthesis   
    Iw = I9                         #2.187 ns
    A = AcropFast(T2, C1, Iw, I9) # include all functions 1.023 ms

    # Control variables this is a temporary implementation no real time 6.467 ns
    U1, U2, U3, U4, U5, U7, U9, U10, U12 =  
     U1f(t), U2f(t), U3f(t), U4f(t), U5f(t),  
     U7f(t), U9f(t), U10f(t), U12f(t)

    # Side and roof control using pBand, ofset
    ofset_t = ofset(t)          #30.078 ns
    pBand_t = VentpBand(t)      #30.078 ns
    U11 = Tset(t)               #30.078 ns

    U6, U8 = vent_control(U11, T2, pBand_t, ofset_t) #24.041 ns

    # Control de la tuberia de calentamiento
    U_HeatPipe, integral_error = pi_control(U11, T2, u[5], 0.075, 0.5, 0.5, 0.1) #24.614 ns
    I3 = U_HeatPipe* HEAT_PIPE + (1 - U_HeatPipe) * T2 # Pipe temperature #300 ns
 
     # compute RHS
    #du[1] = rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12) #u"K * s^-1" 
    #du[2] = rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)#u"K * s^-1" 
    #du[3] = rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10) # u"mg * m^-3 * s^-1"
    #du[4] = rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7,U8, U9)#u"Pa* s^-1"

    du[1], du[2], du[3], du[4] = rhs_fast(T1, T2, C1, V1, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A, U1,  U2, U3,  U4, U5, U6, U7,U8, U9,U10, U12)
    # 5.098 ms

    du[5] = U11 - T2 #32ns
    u[5] = integral_error # 2.424 ns 
    u[6] = I3 # 2.424 ns 
    u[7] = U6 # 2.424 ns 
    u[8] = U8 # 2.424 ns 
    u[9] = A # 2.424 ns 
    #return [dy1, dy2, dy3, dy4]
end



#@everywhere 
function myCallback(sol, t, integrator)
    # Perform additional computations here
end    


# Main simulation loop
@btime begin
    # Days of simulation
    No_days = 90
    days = range(0, stop = No_days, length = No_days+1)
    # Initial conditions
    # y = [T1, T2, C1, V1, integral_error0, I3, U6, U8] 
    
    # Initial date
    p0 = Dates.datetime2unix(start_date)
    
    #Array to store solution 
    global u = [ 18+273.15, 23+273.15, 575 , 1200, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    global time_global = p0

    # Time span: 24 hours in seconds starting form initial time
    t_span = (0, 24 * 3600.0)

    # Define time evaluation points: every 10 minutes (600 seconds)
    global t_eval = (t_span[1]:600:t_span[2])  # Every 600 seconds (10 minutes)

    println("Inicio loop diario")

    # units [u"K", u"K", u"mg * m^-3", u"Pa", none, u"K", u"nounits",u"nounits" ]
    global u0 = [ 18+273.15, 23+273.15, 610 , 1200, 0.0, 273.15, 0, 0, 0 ]
        
    for day_i in days
        # Date as UnixTime for weather data 
        p = p0 + 86400*day_i

        # Define the ODE problem, passing `nothing` as parameters
        prob = ODEProblem(rhsPLL_PID!, u0, t_span, p)

        # Solve the problem with time evaluation points
        #println("Solucion dia $day_i ")
        sol = solve(prob, BS3(), saveat=t_eval)

        # Append solution and time         
        global u = hcat(u, sol[1:9,:])
        global time_global = vcat(time_global, p0.+t_eval .+ day_i*86400)
        global u0 = sol[1:9,end]
    end
end


TT=Tset(vcat(0,t_eval))
TTout=Toutf(time_global)
RH_c = rhf(u[2,:],u[4,:])
I2T = Iglobalf(time_global)

plot((time_global.-time_global[1])/(60*60),[u[1,:].-273.15, u[2,:].-273.15, u[6,:].-273.15, TTout.-273.15], title="Greenhouse T1 & T2", legend=:topright, ylims=(0,100))



time_axis = ((time_global.-time_global[1])/(60*60)).%24



l = @layout [a b; c d]

p1=plot()
p2=plot()
p3=plot()
p4=plot()
plot!(p1,time_axis, u[1,:].-273.15, label="T_{canopy}", lw=1, ylims=(0,55))        # First series with labe
plot!(p1,time_axis, u[2,:].-273.15, label="T_{Air}", lw=1, ylims=(0,55))    # Add second series, dashed line
plot!(p1,time_axis, TTout.-273.15, label="T_{out}", lw=1, ylims=(0,55), alpha=0.7)   # Add third series, dotted line
title!("Temperature")
xlabel!("t (600 s)")
ylabel!("C")
#plot!(time_axis, I2T.*100/800, label="Iglobal")
plot!(p2,time_axis, u[4,:], ylims=(500,2600))
title!("Vapour preasure")
xlabel!("t (600 s)")
ylabel!("Pa")
plot!(p3,time_axis, RH_c)
title!("RH")
ylabel!("%")
plot!(p4,time_axis, u[3,:], ylims=(0,610))
title!(L"CO_2")
ylabel!(L"\frac{mg}{m^3}")
# Customize the plot (optional)



plot(p1, p2, p3, p4, layout = l)
# Show the plot

#p2= plot(time_axis, u[3,:], label="C1_{canopy}", lw=2)        # First series with label

# Show the plot
#display(p2)

#p2 = plot([RH_c, 0*RH_c.+100], ylim=[0,150])
# Show the plot
#display(p2)


#p3= plot(time_axis, u[4,:], label="PV_{canopy}", lw=2)        # First series with label
# Show the plot
#display(p3)

