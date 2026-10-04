using CSV
using DataFrames
using Dates
using Unitful
using DifferentialEquations
using BenchmarkTools



# Include the script that imports constants
include("../utils/import_constants_NoUnits.jl")

# include all functions
include("../utils/climate_functionsV2_NoUnits.jl")

# Import weather data

start_date =  Dates.DateTime(1998, 7, 11, 8, 0)  # Start date and time
end_date   =  Dates.DateTime(1998, 7, 12, 10, 0)  # End date and time

include("../utils/readweatherdata_NoUnits.jl")

include("../utils/temp_controls.jl")

# Initial conditions
# y = [T1, T2, C1, V1] units [u"K", u"K", u"mg * m^-3", u"Pa"]
y0 = [18+273.15, 23+273.15, 680 , 1200]

t = Dates.datetime2unix(start_date)+ 600 # time in seconds

#@btime  dy = rhs!(y0, t)
#dy = rhs!(y0, t)
#print(dy)

println("Inicio solucion de sistema de ODE para Clima del invernadero")

# Time span: 24 hours in seconds starting form initial time
t_span = (0, 24 * 3600.0)

p = Dates.datetime2unix(start_date)
# Define the ODE problem, passing `nothing` as parameters
prob = ODEProblem(rhs!, y0, t_span, p)


# Define time evaluation points: every 10 minutes (600 seconds)
t_eval = (t_span[1]:600:t_span[2])  # Every 600 seconds (10 minutes)


# Measure the time taken to solve the ODE using BenchmarkTools
@btime solve($prob, FBDF(), saveat=$t_eval)
# Trapezoid()  8.571 ms (333725 allocations: 5.19 MiB)
# QNDF()      23.523 ms (819563 allocations: 12.65 MiB)
# BS3()       37.827 ms (1499374 allocations: 23.06 MiB)
# FBDF()      48.172 ms (1689795 allocations: 26.02 MiB)
# Tsit5()     53.546 ms (2105611 allocations: 32.38 MiB)
# DP5()       61.941 ms (2325378 allocations: 35.75 MiB)
# Vern7()     82.261 ms (3332109 allocations: 51.30 MiB)



# Solve the problem with time evaluation points
sol = solve(prob, DP5(),saveat=t_eval)



plot(sol, title="Greenhouse", legend=:topright)
