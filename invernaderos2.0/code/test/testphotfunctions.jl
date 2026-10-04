using CSV
using DataFrames
using Dates
using Unitful
using BenchmarkTools


# Include the script that imports constants
include("../utils/import_constants.jl")

# include all functions
include("../utils/cropphoto_functions.jl")

# Define reazonable values of all input variables 
T1 =  uconvert(u"K", 20u"°C")
T2 = uconvert(u"K", 22u"°C")
C1 = 680u"mg * m^-3"
V1 = 1200u"Pa"
I2 = 1000u"W * m^-2" 

# Test A_B
#Set fixed value for I
I = 0.5*J_max
C_values = range(1u"mg*m^-3", stop=2500u"mg*m^-3", length=50)
using Plots
plot(title="Function A_B vs C for several T", xlabel="C", ylabel="A_B(T, C, I)", lw=2)
#  Loop over each value of T and plot the corresponding function
A_B_values1 = [A_B(T1, C, 1000e-6u"mol * m^-2 * s^-1") for C in C_values]
A_B_values2 = [A_B(T2, C, 500e-6u"mol * m^-2 * s^-1") for C in C_values]
plot(C_values, [A_B_values1, A_B_values2]/1e-6)  # Plot each T with a different label

println(Ar_P(T1))
println(Ar_C(T1,C1))

Ar= Aresult(T1,C1,J_max*10)
println(Ar)

Ac = Acrop(T1, C1, I2, 2)
println(Ac)