using CSV
using DataFrames
using Dates
using Unitful
using BenchmarkTools


# Include the script that imports constants
include("../utils/import_constants_NoUnits.jl")

# include all functions
include("../utils/cropphoto_functions_NoUnits.jl")

# Define reazonable values of all input variables 
T1 =  20 + 273.15 #uconvert(u"K", 20u"°C")
T2 = 22 + 273.15 #uconvert(u"K", 22u"°C")
C1 = 680#u"mg * m^-3"
V1 = 1200#u"Pa"
I2 = 1000#u"W * m^-2" 


# Test A_B
#Set fixed value for I
I = 0.5*J_max
C_values = range(1, stop=2500, length=50)
using Plots
plot(title="Function A_B vs C for several T", xlabel="C", ylabel="A_B(T, C, I)", lw=2)
#  Loop over each value of T and plot the corresponding function
A_B_values1 = [A_B(T1, C, 1000e-6 ) for C in C_values]
A_B_values2 = [A_B(T2, C, 500e-6 ) for C in C_values]
plot(C_values, [A_B_values1, A_B_values2]/1e-6)  # Plot each T with a different label



println(Ar_P(T1))
println(Ar_C(T1,C1))

Ar= Aresult(T1, C1, J_max*10)
println(Ar)

println("Optimizacion Ar_P")
@btime (1.5*V_cmax(T1)-Rd_day)
@btime V_cmax(T1)
@btime Ar_P(T1)

println("Optimizacion Ar_C")
@btime Ar_C(T1,C1)
println("Optimizacion Ar_J")
@btime Ar_J(T1,C1,I2)

println("Test A crop normal")
@btime Ac = Acrop(T1, C1, I2, 2)
#println(Ac)

println("Test A crop fast")
@btime Ac = AcropFast(T1, C1, I2, 2)