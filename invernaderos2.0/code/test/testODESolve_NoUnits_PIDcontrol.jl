using CSV
using DataFrames
using Dates
using Unitful
using DifferentialEquations
using Distributed
using BenchmarkTools



# Include the script that imports constants
include("../utils/import_constants_NoUnits.jl")

# include all functions
include("../utils/climate_functionsV2_NoUnits.jl")

# Import weather data

start_date =  Dates.DateTime(1998, 7, 11, 8, 0)  # Start date and time
end_date   =  Dates.DateTime(1998, 11, 12, 10, 0)  # End date and time

include("../utils/readweatherdata_NoUnits.jl")

include("../utils/temp_controls.jl")

include("../utils/control_functions.jl")

# Define the right-hand side function of the ODE
@everywhere function rhsPLL_PID!(du, u, p, t)
    # p[1] = start_date(unix)
   

    # Unpack  y = [T1, T2, C1, V1] and u[5] = integral_error 
    # unpack y 
    T1, T2, C1, V1 = u[1], u[2], u[3], u[4]   
    integral_error = u[5]

    # from simulation time t to date(unix) time t_calendar
    t_calendar = t  + p

    #weather data is given in calendar days 
    I1, I2, I4, I5, I6, I7, I8, I9, I10, I11 = 
     LAIf(t_calendar), Iglobalf(t_calendar),  
     Tskyf(t_calendar), Toutf(t_calendar), TmechCoolf(t_calendar), 
     Tsoilf(t_calendar), WindSpeedf(t_calendar), Idocelf(t_calendar),  
     CO2outf(t_calendar),  VPoutf(t_calendar)

    #Assimilates function photosynthesis   
    A = Af(t) 

    # Control variables this is a temporary implementation no real time
    U1, U2, U3, U4, U5, U7, U9, U10, U11, U12 =  
     U1f(t), U2f(t), U3f(t), U4f(t), U5f(t),  
     U7f(t), U9f(t), U10f(t), U11f(t), U12f(t)

    # Side and roof control using pBand, ofset
    U6, U8 = vent_control(U11, T2, 15, 1)

    # Control de la tuberia de calentamiento
    U_HeatPipe, integral_error = pi_control(U11, T2, u[5], 0.075, 0.5, 0.029, 0.1)
    I3 = U_HeatPipe * HEAT_PIPE + (1 - U_HeatPipe) * T2 # Pipe temperature 
   
    

 
     # compute RHS
    du[1] = rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12) #u"K * s^-1" 
    du[2] = rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)#u"K * s^-1" 
    du[3] = rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10) # u"mg * m^-3 * s^-1"
    du[4] = rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7,U8, U9)#u"Pa* s^-1"

    du[5] = U11 - T2
    u[5] = integral_error 
    u[6] = I3
    u[7] = U6
    u[8] = U8
    #return [dy1, dy2, dy3, dy4]
end



@everywhere function myCallback(sol, t, integrator)
    # Perform additional computations here
end    


@everywhere begin
# Initial conditions
# y = [T1, T2, C1, V1] units [u"K", u"K", u"mg * m^-3", u"Pa"]
    y0 = [18+273.15, 23+273.15, 575 , 1200]
    integral_error0 = 0.0  # Initial integral of the error
    I3r = 275.15

    u0 = [ y0[1], y0[2], y0[3], y0[4], integral_error0, I3r, 0.0, 0.0 ]
    t = Dates.datetime2unix(start_date) + 600 # time in seconds

    println("Inicio solucion de sistema de ODE para Clima del invernadero")

# Time span: 24 hours in seconds starting form initial time
    t_span = (0, 24 * 3600.0*60)
    p = Dates.datetime2unix(start_date)
# Define the ODE problem, passing `nothing` as parameters
    prob = ODEProblem(rhsPLL_PID!, u0, t_span, p)


# Define time evaluation points: every 10 minutes (600 seconds)
    t_eval = (t_span[1]:300:t_span[2])  # Every 600 seconds (10 minutes)

#@btime sol = solve(prob, BS3(),saveat=t_eval)
# Solve the problem with time evaluation points
    sol = solve(prob,Tsit5() , saveat=t_eval)
end


plot(sol[1,:].-273.15, title="Greenhouse T1", legend=:topright, ylims=(20,40))
#plot(sol[2,:].-273.15, title="Greenhouse T2", legend=:topright)
#plot(sol[3,:], title="Greenhouse C1", legend=:topright)
#plot(sol[4,:], title="Greenhouse VP", legend=:topright)
plot(sol[6,:].-273.15, title="Greenhouse I3", legend=:topright, ylims=(20,100))

println(minimum(sol[2,:])-273.15)
println(mean(sol[2,:])-273.15)
println(maximum(sol[2,:])-273.15)
plot(t_eval, [sol[1,:].-273.15,sol[2,:].-273.15,sol[6,:].-273.15], 
    title="Greenhouse T1", legend=:topright, ylims=(10,100))