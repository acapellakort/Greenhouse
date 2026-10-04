# In this version rhsPLL_PID! is updated to save all control values
# Probably this code will generate all data and save it as experimental 
# measurements. Th inference is made elsewhere. 


# Required packages
using CSV
using DataFrames
using Unitful
using Dates
using LaTeXStrings

using SparseArrays
using Statistics
using DifferentialEquations
#using Distributed
#using BenchmarkTools
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

# Setting of start and end date-time
# Initial hour 0 so that the time stamp is sincronized with the controls
start_date  =  Dates.DateTime(2018, 8, 15, 0, 0)                         # Start date  
start_time = Dates.datetime2unix(start_date)                             # Star time (same as start date but in UNIX time)
No_days     = 2                                                          # Number of days to simulate
end_date    =  Dates.DateTime(2018, 8, 15 + No_days, 10, 0)              # End date and time
#end_date   =  Dates.DateTime(2018, 9, 15, 10, 0)                        # End date and time

# Load experimental data
# Import weather data between start_date and end_date for the experiment
include("./utils/readweatherdata_experimento.jl")

# Load experiment observation to a csv-file to a dataframe 
datafilename = "real_observed_data.csv"
observed_data = DataFrame(CSV.File(datafilename))

# Filter experimental data for missing values 
observed_data = filter(row -> !ismissing(row.time) && row.time > Dates.datetime2unix(start_date) && row.time < Dates.datetime2unix(end_date), observed_data)

# Observation for the logLikelihood
tcan_data   = observed_data.T1 
tair_data   = observed_data.T2 
rh_data     = observed_data.RH  
co2_data    = observed_data.C1 
time_data   = observed_data.time 
# Number of data registries (rows in observed_data after filtering for misssing values)
n = length(observed_data.T1) 

# Assumed standard deviation error for the observations (taken in consult aaron) usded in the logLikelihood
error_std   = 0.5*[0.5, 0.5, 5, 0.1,  0.1 ] 

##########################################################################

# plot inference results
function graficas(mcmc, QoIval, k, savefigs = true)
    # Plot for -energy of the mcmc
    p0=plot()
    plot!(p0, mcmc.Output[:,end])
    
    #plot of the variable k
    p=plot()
    histogram!(p,  mcmc.Output[burn_in:end,k], color=:blue, bins=10)
    #vline!(p, [QoI_dict[QoIval][2]], lw = 3, linecolor = :green)

    xlabel!(QoIval)
    plot!(p,legend = :false)

    if savefigs == true
        png(p0, "figures/"*QoIval*"_chain")
        png(p, "figures/"*QoIval*"_posterior")
    end
end

# Control functions
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
    I1 = LAIf(t_calendar)*0.1
    I2 = Iglobalf(t_calendar)
    I4 = Tskyf(t_calendar) 
    I5 = Toutf(t_calendar) +2.0
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

# parameter for  inference
QoI_dict = Dict(
    "alpha1"    => [1000.0, 1500.0, 2000.0 ],  
    "alpha2"    => [   0.1,   0.35,    0.99],
    "alpha4"    => [   1.0,    5.0,   50.0 ],
    "beta2"     => [   0.1,    0.7,    2.0 ],
    "gamma3"    => [ 100.0,  275.0,  500.0 ],
    "gamma4"    => [  20.0,   82.0,  200.0 ],
    "nu4"       => [  1e-5,   1e-4,   1.0  ],
    "nu4CO2"    => [  1e-7,   1e-3,  1e-4  ],
    "nu7"       => [   0.1,   0.85,    2.0 ],
    "eta1"      => [  0.05,   0.10,   0.95 ],
    "tau1"      => [  0.05,   0.77,   0.99 ], 
    "psi2"      => [  50.0, 13300.0, 25000 ]
    ) 


# Inferrence setup
global QoI = ["alpha1", "alpha2", "alpha4",  "beta2", "gamma3", "gamma4", "nu4", "nu4CO2", "nu7", "eta1", "tau1", "psi2" ]   
chain_length = 3000
burn_in      = 1500

# Traking print
println("Parameter(s) to infer : ")
println(QoI)

# Initial conditions for the ODE system
vp0 = rh_data[1]/100*611.2*exp(  17.67*(tair_data[1]-273.15)/(tair_data[1]-29.65))

global u0s = [tair_data[1], tair_data[1], co2_data[1], vp0, 0.0, 0.0]    # units [u"K", u"K", u"mg * m^-3", u"Pa", u"K", u"K * s" ]             

# 
global U1s, U2s, U3s, U4s, U5s, U6s, U7s, U8s, U9s, U10s, U11s, U12s, I3s = control4inference(datafilename)

# Forward map 
function FM_SimLoop!(x, theta_name)
# Forward Map simulation 

    # Set the QoI value in forward simulations
    n = length(theta_name)
    for k in 1:n
        var_name = theta_name[k]
        var_value = x[k]
        global_var_name = Symbol(var_name)
        eval(:(global $global_var_name = $var_value))
    end 

    # Solve de ODE for the full period of time
    p = start_time
    # Experimento define values to compare
    t_eval = time_data.-time_data[1]
    t_span = (t_eval[1],t_eval[end])
    # Define the ODE problem with current parameters `p`
    prob = ODEProblem(rhsPLL_PID_sim!, u0s, t_span, p) 
    # Solve the ODE with time evaluation points
    sol = solve(prob, BS3(), saveat=t_eval)
    # Report solution
    RH_obs      = rhf(sol[2, :], sol[4, :])
    CO2_obs     = sol[3,:]
    T1_canopy   = sol[1,:]
    T2_air      = sol[2,:]

    return  T1_canopy, T2_air, RH_obs, CO2_obs
end

# Precition matrix 
function create_tridiagonal_matrix(n, v)
    # Create a sparse matrix of size n x n
    A = spzeros(n, n)

    # Fill in the main diagonal, first off-diagonals, and second off-diagonals
    for i in 1:n
        A[i, i] = v[1]           # Main diagonal
        if i > 1
            A[i, i - 1] = v[2]  # First off-diagonal (sub-diagonal)
            A[i - 1, i] = v[2]  # First off-diagonal (super-diagonal)
        end
        if i > 2
            A[i, i - 2] = v[3]   # Second off-diagonal (sub-diagonal)
            A[i - 2, i] = v[3]   # Second off-diagonal (super-diagonal)
        end
    end

    return A
end

function PriorSupp(x)
# Prior's support    
    n = length(x)
    X = true
    for k in 1:n
        X = X && (QoI_dict[QoI[k]][1] .< x[k] &&  x[k]< QoI_dict[QoI[k]][3])
    end
    return X
end

# Correlation matrices for loklikelihood
n_obs = length(tair_data)
SigmaT = create_tridiagonal_matrix(n_obs,[1, -0.5, -0.026] )
SigmaRH = create_tridiagonal_matrix(n_obs, [1,0.7, 0.25])
SigmaCO2 = create_tridiagonal_matrix(n_obs, [1, 0.03, 0.01])

# energy (uniform prior)
function energy(x)
    tcan_p, tair_p, rh_p, co2_p= FM_SimLoop!(x, QoI)
    logLikelihood = #- sum((tcan_data.-tcan_p).^2)./(error_std[1])^2
                    #- sum((tair_data.-tair_p).^2)./(error_std[2]).^2
                    #- sum((rh_data .- rh_p).^2  )./(error_std[3]).^2
                    #- sum((co2_data.-co2_p).^2  )./(error_std[4])^2
                    -transpose((tair_data.-tair_p)) * (SigmaT * (tair_data.-tair_p))/(error_std[2]).^2
                    -transpose((rh_data .- rh_p)) * (SigmaRH * (rh_data .- rh_p))/(error_std[3]).^2
                    -transpose((co2_data.-co2_p)) * (SigmaCO2 * (co2_data.-co2_p))/(error_std[4])^2
                    #-((tair_data.-tair_p)\SigmaT)*(tair_data.-tair_p)/(error_std[2]).^2
                    #-((rh_data .- rh_p)\SigmaRH)*(rh_data .- rh_p)/(error_std[3]).^2
                    #-((co2_data.-co2_p)\SigmaCO2)*(co2_data.-co2_p)/(error_std[4])^2
                    +0.85*sum((tair_data.-tair_p).*(rh_data .- rh_p))./(error_std[2]*error_std[3])
                    +0.1*sum((tair_data.-tair_p).*(co2_data.-co2_p))./(error_std[2]*error_std[4])
                    +0.27*sum((rh_data .- rh_p).*(co2_data.-co2_p))./(error_std[3]*error_std[4])

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



# Sample data (replace with your actual data)
data = DataFrame(Temperature=tair_data,
                 RH=rh_data,
                 CO2=co2_data)

# mcmc run

correlation_matrix = cor(Matrix(data))
println(correlation_matrix)

println("Inicio MCMC")
@time Run!(QoI_mcmc, T=chain_length, x0=x0 , xp0=xp0 )

# plot results


for i in 1:n
    graficas(QoI_mcmc, QoI[i], i, true)
end

loglike, indx = findmin(QoI_mcmc.Output[:,end])
map = QoI_mcmc.Output[indx,1:end-1]
postmean = mean.( eachcol(QoI_mcmc.Output[burn_in:end,1:end-1]))
v_map = FM_SimLoop!(map,QoI)
v_mean = FM_SimLoop!(postmean,QoI)

plot(v_map[2].-273.15)
plot!(v_mean[2].-273.15)
plot!(tair_data.-273.15)

plot(v_map[3])
plot!(v_mean[3])
plot!(rh_data)

plot(v_map[4])
plot!(v_mean[4])
plot!(co2_data)

xt = ones(n)
for k in 1:n
    xt[k]= QoI_dict[QoI[k]][2]
end

vm = FM_SimLoop!(xt,QoI) 
plot(vm[4])