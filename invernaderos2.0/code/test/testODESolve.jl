using CSV
using DataFrames
using Dates
using Unitful
using DifferentialEquations
using BenchmarkTools



# Include the script that imports constants
include("../utils/import_constants.jl")

# include all functions
include("../utils/climate_functionsV2.jl")

# Import weather data

start_date =  Dates.DateTime(1998, 7, 11, 8, 0)  # Start date and time
end_date   =  Dates.DateTime(1998, 7, 12, 10, 0)  # End date and time

include("../utils/readweatherdata.jl")

include("../utils/temp_controls.jl")

# Initial conditions
# y = [T1, T2, C1, V1]
y0 = [uconvert(u"K " ,18u"°C"), uconvert(u"K " ,23u"°C"), 680u"mg * m^-3", 1200u"Pa"]

t = Dates.datetime2unix(start_date)*u"s" + 600u"s"

#@btime  dy = rhs!(y0, t)
#dy = rhs!(y0, t)
#print(dy)

println("Inicio solucion de sistema de ODE para Clima del invernadero")

# Time span: 24 hours in seconds starting form initial time
t_span = (0u"s", 24 * 3600.0u"s")

p = Dates.datetime2unix(start_date)u"s"
# Define the ODE problem, passing `nothing` as parameters
prob = ODEProblem(rhs!, y0, t_span, p)


# Define time evaluation points: every 10 minutes (600 seconds)
t_eval = (t_span[1]:300u"s":t_span[2])  # Every 600 seconds (10 minutes)

# Solve the problem with time evaluation points
sol = solve(prob,saveat=t_eval)


plot(sol, title="Greenhouse", legend=:topright)
