# In this version rhsPLL_PID! is updated to save all control values
# Probably this code will generate all data and save it as experimental 
# measurements. Th inference is made elsewhere. 


# Required packages
using CSV
using DataFrames
using Dates
using Unitful
using DifferentialEquations
using OrdinaryDiffEq
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

# Main simulation loop 
# Initial hour 0 so that the time stamp is sincronized with the controls
start_date =  Dates.DateTime(1998, 7, 11, 0, 0)  # Start date and time
end_date   =  Dates.DateTime(1998, 11, 12, 10, 0)  # End date and time

# ############# Eventually the dates should come from the experimental data<--------

# Import weather data between start_date and end_date
include("./utils/readweatherdata_NoUnits.jl")

No_days = 2-1  # Number of days to simulate
days = range(0, stop = No_days, length = No_days + 1)
start_time = Dates.datetime2unix(start_date)

# Time span for each day in seconds (24 hours)
t_span = (0, 24 * 3600.0)

# Define time evaluation points every 10 minutes (600 seconds)
global rep_time = 3600 # reporting time  
global t_eval = (t_span[1]:rep_time:t_span[2])

############# Load experiment csv-file to a dataframe 
datafilename = "../observed_data.csv"
observed_data = DataFrame(CSV.File(datafilename))

# Lenght of the observed data
n = length(observed_data.T1) 

# add error to the observations 
error_std   = 1.5*[0.5, 0.5, 5, 0.1,  0.1 ] 
tcan_data   = observed_data.T1 
tair_data   = observed_data.T2 
rh_data     = observed_data.RH  
co2_data    = observed_data.C1 

##########################################################################

# plot inference results
function graficas(mcmc, QoIval, k, savefigs = true)
    p0=plot()
    plot!(p0, mcmc.Output[:,end])
    
    p=plot()
    histogram!(p,  mcmc.Output[burn_in:end,k], color=:blue, bins=10)
    vline!(p, [QoI_dict[QoIval][2]], lw = 3, linecolor = :green)

    xlabel!(QoIval)
    plot!(p,legend = :false)

    if savefigs == true
        png(p0, "figures/"*QoIval*"_chain")
        png(p, "figures/"*QoIval*"_posterior")
    end
end

# define control functions
function control4inference(filename)
    # This function read and interpolate the controls used in the measurements
    # All the measured data is contained in the file = filename 
    # Data is organized as follows:
    #
    # Time, T1, T2, C1, V1, RH, Photo, Tpipe, U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12  
    #
    df = DataFrame( CSV.File(filename))
    time_abs = df.time.-df.time[1]
    U1d = DataInterpolations.ConstantInterpolation(df.U1,time_abs )
    U2d = DataInterpolations.ConstantInterpolation(df.U2,time_abs )
    U3d = DataInterpolations.ConstantInterpolation(df.U3,time_abs )
    U4d = DataInterpolations.ConstantInterpolation(df.U4,time_abs )
    U5d = DataInterpolations.ConstantInterpolation(df.U5,time_abs )
    U6d = DataInterpolations.ConstantInterpolation(df.U6,time_abs )
    U7d = DataInterpolations.ConstantInterpolation(df.U7,time_abs )
    U8d = DataInterpolations.ConstantInterpolation(df.U8,time_abs )
    U9d = DataInterpolations.ConstantInterpolation(df.U9,time_abs )
    U10d = DataInterpolations.ConstantInterpolation(df.U10,time_abs )
    U11d = DataInterpolations.ConstantInterpolation(df.U11,time_abs )
    U12d = DataInterpolations.ConstantInterpolation(df.U12,time_abs )
    I3d = DataInterpolations.ConstantInterpolation(df.Tpipe,time_abs )
    return U1d, U2d, U3d, U4d, U5d, U6d, U7d, U8d, U9d, U10d, U11d, U12d, I3d
end

# Define the right-hand for the inference loop 
function rhsPLL_PID_sim!(du, u, p, t)
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
    U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12 = U1s(t), U2s(t), U3s(t), U4s(t), U5s(t), U6s(t), U7s(t), U8s(t), U9s(t), U10s(t), U11s(t), U12s(t)
    
    # Heat pipe temperature
    I3 = I3s(t)  

    # Falta control de la fuente de CO2

    # Compute right-hand sides (RHS) of the ODEs
    du[1], du[2], du[3], du[4] = rhs_fast(T1, T2, C1, V1, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A, U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12)
    
    # Update integral error and additional variables
    du[5] = U11 - T2        # Integral of error for PI control
    u[5]  = integral_error   
    u[6]  = A                # Photosynthesis response
    
end

##########################################################################
#                                                                        #
#       Set up of the inference                                          #
#                                                                        #
##########################################################################

# parameter for  inference
QoI_dict = Dict(
    "alpha1"    => [1000.0, 1500.0, 2000.0 ],  
    "alpha2"    => [   0.1,   0.35,    0.99],
    "alpha4"    => [   1.0,    5.0,   50.0 ],
    "beta2"     => [   0.1,    0.7,    2.0 ],
    "gamma3"    => [ 100.0,  275.0,  500.0 ],
    "gamma4"    => [  20.0,   82.0,  200.0 ],
    "nu4"       => [  1e-5,   1e-2,   1e-1 ],
    "nu7"       => [   0.1,   0.85,    2.0 ],
    "eta1"      => [  0.05,   0.10,   0.95 ],
    "tau1"      => [  0.05,   0.77,   0.99 ]
    ) 

# Inferrence setup
global QoI = ["beta2", "gamma3", "gamma4", "nu4"]   
chain_length = 25000
burn_in      = 20000


println("Parameter(s) to infer : ")
println(QoI)

# Initial conditions for the ODE system
global u0s = [18+273.15, 23+273.15, 575, 1200, 0.0, 0.0]
# units [u"K", u"K", u"mg * m^-3", u"Pa", u"K", u"K * s" ]
start_time = Dates.datetime2unix(start_date)
# Arrays to store the solution and time variable
global us = u0s
global time_globals = start_time
# Time span for each day in seconds (24 hours)
t_span = (0, 24 * 3600.0)
# 
global U1s, U2s, U3s, U4s, U5s, U6s, U7s, U8s, U9s, U10s, U11s, U12s, I3s = control4inference(datafilename)

# Forward map 
function FM_SimLoop!(x, theta_name)

    # Initial condition
    us = u0s
    # Set the QoI value in forward simulations
    n = length(theta_name)
    for k in 1:n
        var_name = theta_name[k]
        var_value = x[k]
        global_var_name = Symbol(var_name)
        eval(:(global $global_var_name = $var_value))
    end 

    # Day's loop
    for day_i in days
        p = start_time + 86400 * day_i  # Convert day index to Unix time
        
        # Define the ODE problem with current parameters `p`
        prob = ODEProblem(rhsPLL_PID_sim!, u0s, t_span, p)
        
        # Solve the ODE with time evaluation points
        sol = solve(prob, Tsit5(), saveat=t_eval)
        
        # Append the solution and update initial conditions for next day
        us = hcat(us, sol[1:6,:])

        #global time_global = vcat(time_global, start_time.+t_eval .+ day_i*86400)
        global u0s = sol[1:6,end]
    end

    RH_obs      = rhf(us[2, :], us[4, :])
    CO2_obs     = us[3,:]
    T1_canopy   = us[1,:]
    T2_air      = us[2,:]
    return  T1_canopy, T2_air, RH_obs, CO2_obs
end

# Prior support
function PriorSupp(x)
    n = length(x)
    X = true
    for k in 1:n
        X = X && (QoI_dict[QoI[k]][1] .< x[k] &&  x[k]< QoI_dict[QoI[k]][3])
    end
    return X
end

# energy (uniform prior)
function energy(x)
    tcan_p, tair_p, rh_p, co2_p= FM_SimLoop!(x, QoI)
    logLikelihood = - sum((tcan_data.-tcan_p).^2)./(error_std[1])^2
                    - sum((tair_data.-tair_p).^2)./(error_std[2]).^2
                    - sum((rh_data .- rh_p).^2  )./(error_std[3]).^2
                    - sum((co2_data.-co2_p).^2  )./(error_std[4])^2
    return -logLikelihood
end

n = length(QoI)
#Define mcmc object 
QoI_mcmc = jtwalk( n=n, U=energy, Supp=PriorSupp)

# Draw random Initial points from uniform prior
x0 = ones(n)
xp0 = ones(n)
for k in 1:n
    x0[k], xp0[k] = rand(2).*(QoI_dict[QoI[k]][3].-QoI_dict[QoI[k]][1]).+QoI_dict[QoI[k]][1]
end

# mcmc run

println("Inicio MCMC")
Run!(QoI_mcmc, T=chain_length, x0=x0 , xp0=xp0 )

# plot and 
graficas(QoI_mcmc, QoI[1], 1, true)
graficas(QoI_mcmc, QoI[2], 2, true)