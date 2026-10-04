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
end_date   =  Dates.DateTime(1998, 10, 20, 10, 0)  # End date and time

include("../utils/readweatherdata_NoUnits.jl")

include("../utils/temp_controls.jl")




# Define the right-hand side function of the ODE
@everywhere function rhsPLL!(dy, y, p, t)
    # In this function y = [T1, T2, C1, V1])
        #Tout = I5
        #WindSpeed = I8
        #CO2out = I10 
        #Iglobal = I2   
        #VPout = I11
        # I1 = LAI  
        # I2 = Iglobal  
        # I3 = Tpipef    
        # I4 = Tskif
        # I5 = Tout
        # I6 = TmechCool  
        # I7 = Tsoil 
        # I8 = WindSpeed m/s
        # I9 = Idocel
        #I10 = CO2out  
        #I11 = VPout

    
    # Move evaluation of weather to current date
    t_floatR = t  + p
    #weather data
    I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11 = 
     LAIf(t_floatR), Iglobalf(t_floatR),  Tpipef(t_floatR), 
     Tskyf(t_floatR), Toutf(t_floatR), TmechCoolf(t_floatR), 
     Tsoilf(t_floatR), WindSpeedf(t_floatR), Idocelf(t_floatR),  
     CO2outf(t_floatR),  VPoutf(t_floatR)

    # Control variables this is a temporary implementation no real time
    U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12 =  
    U1f(t),  U2f(t), U3f(t), U4f(t), U5f(t), U6f(t), 
    U7f(t), U8f(t), U9f(t), U10f(t), U12f(t)
    
    #Assimilates function photosynthesis   
    A = Af(t)    
    
    # unpack y 
    T1, T2, C1, V1 = y[1], y[2], y[3], y[4]

     # compute RHS
    dy[1] = rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12) #u"K * s^-1" 
    dy[2] = rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)#u"K * s^-1" 
    dy[3] = rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10) # u"mg * m^-3 * s^-1"
    dy[4] = rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7,U8, U9)#u"Pa* s^-1"

    #return [dy1, dy2, dy3, dy4]
end



@everywhere function myCallback(sol, t, integrator)
    # Perform additional computations here
end    


@everywhere begin
# Initial conditions
# y = [T1, T2, C1, V1] units [u"K", u"K", u"mg * m^-3", u"Pa"]
    y0 = [18+273.15, 23+273.15, 575 , 1200]
    t = Dates.datetime2unix(start_date)+ 600 # time in seconds

    println("Inicio solucion de sistema de ODE para Clima del invernadero")

# Time span: 24 hours in seconds starting form initial time
    t_span = (0, 24 * 3600.0*90)
    p = Dates.datetime2unix(start_date)
# Define the ODE problem, passing `nothing` as parameters
    prob = ODEProblem(rhsPLL!, y0, t_span, p)


# Define time evaluation points: every 10 minutes (600 seconds)
    t_eval = (t_span[1]:300:t_span[2])  # Every 600 seconds (10 minutes)


# Measure the time taken to solve the ODE using BenchmarkTools
#    @btime solve($prob, Trapezoid(), saveat=$t_eval)
# Parallel solver trapezoid 8.393 ms (333725 allocations: 5.19 MiB)
# Trapezoid()  8.571 ms (333725 allocations: 5.19 MiB)
# QNDF()      23.523 ms (819563 allocations: 12.65 MiB)
# BS3()       37.827 ms (1499374 allocations: 23.06 MiB)
# FBDF()      48.172 ms (1689795 allocations: 26.02 MiB)
# Tsit5()     53.546 ms (2105611 allocations: 32.38 MiB)
# DP5()       61.941 ms (2325378 allocations: 35.75 MiB)
# Vern7()     82.261 ms (3332109 allocations: 51.30 MiB)



# Solve the problem with time evaluation points
    sol = solve(prob, DP5(),saveat=t_eval)
end


plot(sol[1,:].-273.15, title="Greenhouse T1", legend=:topright)
plot(sol[2,:].-273.15, title="Greenhouse T2", legend=:topright)
plot(sol[3,:], title="Greenhouse C1", legend=:topright)
#plot(sol[4,:], title="Greenhouse VP", legend=:topright)
