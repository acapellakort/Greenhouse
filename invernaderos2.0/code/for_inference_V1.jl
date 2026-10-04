# Required packages
using CSV
using DataFrames
using Dates
using Unitful
using DifferentialEquations
using Distributed
using BenchmarkTools
using LaTeXStrings
using JTwalk

# Import constants used in the simulation
include("./utils/import_constants_NoUnits.jl")

# Include all necessary functions for climate and crop simulations
include("./utils/climate_functionsV3_NoUnits.jl")
# cropphoto_functionsV3_NoUnits depends on climate_functionsV3_NoUnits.jl
include("./utils/cropphoto_functionsV3_NoUnits.jl")

# Import temporary controls (TO BE UPDATED)
include("./utils/temp_controlsV3.jl")

# Include control functions used for the system's response
include("./utils/control_functions.jl")

# Load predefined control settings
include("./utils/import_control_set.jl")

# Define the right-hand side function of the ODE system governing the greenhouse control
function rhsPLL_PID!(du, u, p, t)
    # Unpack variables: [T1, T2, C1, V1, integral_error]
    T1 = u[1]               # Canopy temperature
    T2 = u[2]               # Air temperature
    C1 = u[3]               # CO2 concentration
    V1 = u[4]               # Vapor pressure
    integral_error = u[5]   # Integral error for the PI controller
    
    # Convert simulation time `t` to calendar time
    t_calendar = t + p

    # Retrieve weather and environmental data based on the time
    I1 = LAIf(t_calendar)
    I2 = Iglobalf(t_calendar)
    I4 = Tskyf(t_calendar)
    I5 = Toutf(t_calendar)
    I6 = TmechCoolf(t_calendar)
    I7 = Tsoilf(t_calendar)
    I8 = WindSpeedf(t_calendar)
    I9 = (1 - eta1) * tau1 * eta2 * Idocelf(t_calendar)
    I10 = CO2outf(t_calendar)
    I11 = VPoutf(t_calendar)

    # Photosynthesis assimilation calculation
    Iw = I9
    A = AcropFast(T2, C1, Iw, I9)  # Photosynthesis response

    # Control variables: U1 - U12
    U1, U2, U3, U4, U5, U7, U9, U10, U12 = U1f(t), U2f(t), U3f(t), U4f(t), U5f(t), U7f(t), U9f(t), U10f(t), U12f(t)
    
    # Thermal screen control 
    U1  = (I5 < ToutMax(t) ? 1 : 0) 

    # Heating pipe set point 
    U11 = Tset(t)                        
    
    # Ventilation control using side/roof settings
    U6, U8 = vent_control(U11, T2, VentpBand(t), ofset(t))

    # Forced ventilation control
    U7 = Vent_Forced(U8, U11, T2)

    # PI controller for the heating pipe
    U_HeatPipe, integral_error = pi_control(U11, T2, u[5], 0.075, 0.5, 0.5, 0.1)
    # Heat pipe temperature
    I3 = U_HeatPipe * HEAT_PIPE + (1 - U_HeatPipe) * T2  

    # Light control (on/off)
    U12 = Light_on(t) 

    # Falta control de la fuente de CO2

    # Compute right-hand sides (RHS) of the ODEs
    du[1], du[2], du[3], du[4] = rhs_fast(T1, T2, C1, V1, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A, U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12)
    
    # Update integral error and additional variables
    du[5] = U11 - T2        # Integral of error for PI control
    u[5] = integral_error   
    u[6] = I3               # Heat pipe temperature
    u[7] = U6               # Ventilation control
    u[8] = U8               # Another ventilation control
    u[9] = A                # Photosynthesis response
end

# Main simulation loop 
# Initial hour 0 so that the time stamp is sincronized with the controls
start_date =  Dates.DateTime(1998, 7, 11, 0, 0)  # Start date and time
end_date   =  Dates.DateTime(1998, 11, 12, 10, 0)  # End date and time
# Import weather data between start_date and end_date
include("./utils/readweatherdata_NoUnits.jl")


No_days = 7-1  # Number of days to simulate
days = range(0, stop = No_days, length = No_days + 1)


# Initial conditions for the ODE system
global u0 = [18+273.15, 23+273.15, 575, 1200, 0.0, 0.0, 0.0, 0.0, 0.0]
# units [u"K", u"K", u"mg * m^-3", u"Pa", u"1", u"K", u"1", u"1" ]
start_time = Dates.datetime2unix(start_date)

# Arrays to store the solution and time variable
global u = u0
global time_global = start_time

# Time span for each day in seconds (24 hours)
t_span = (0, 24 * 3600.0)

# Define time evaluation points every 10 minutes (600 seconds)
global t_eval = (t_span[1]:3600:t_span[2])

# De¡finition of forward map
function forwardMap!(x)

    u = u0
    global alpha4 = x
    for day_i in days
        p = start_time + 86400 * day_i  # Convert day index to Unix time
        
        # Define the ODE problem with current parameters `p`
        prob = ODEProblem(rhsPLL_PID!, u0, t_span, p)
        
        # Solve the ODE with time evaluation points
        sol = solve(prob, BS3(), saveat=t_eval)
        
        # Append the solution and update initial conditions for next day
        u = hcat(u, sol[1:9,:])
        #global time_global = vcat(time_global, start_time.+t_eval .+ day_i*86400)
        global u0 = sol[1:9,end]
    end

    RH_obs      = rhf(u[2, :], u[4, :])
    CO2_obs     = u[3,:]
    T1_canopy   = u[1,:]
    T2_air      = u[2,:]
    return  T1_canopy, T2_air, RH_obs, CO2_obs
end



# Observational data simulation

tcan_o, tair_o, rh_o, co2_o= forwardMap!(8)
plot(rh_o)

n = length(tcan_o) 
# Tcan, Tair, RH, CO2
error_std =[0.5, 0.5, 5, 0.1,  0.1 ] 
rh_data = rh_o + error_std[3]*randn(n)
tair_data = tair_o + error_std[1]*randn(n)
tcan_data = tcan_o + error_std[2]*randn(n)

error_std2 =[0.5, 0.5, 0.1, 0.1] 

function energy(x)
    tcan_p, tair_p, rh_p, co2_p= forwardMap!(x[1])
    logLikelihood = -sum((rh_data.-rh_p).^2)./(n*error_std[3]).^2-sum((tcan_data.-tcan_p).^2)./(n*error_std[1])^2
    return -logLikelihood
end

function PriorSupp(x)
    if 1 .< x[1] &&  x[1]< 30 
        return true
    else
        return false
    end
end


ConvectHeatCanopy = jtwalk( n=1, U=energy, Supp=PriorSupp)

println("Inicio MCMC")
Run!(ConvectHeatCanopy, T=5000, x0=6*ones(1), xp0=11.0*ones(1))

plot(ConvectHeatCanopy.Output[:,end])