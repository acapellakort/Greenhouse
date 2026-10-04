
using JSON
using DifferentialEquations
using BenchmarkTools
using Plots

# Define the path to the JSON file
json_file_path = "/Users/usuario/Documents/trabajo/academia/unam/proyectos/invernaderos/invernaderos2.0/code/configfiles/constants-test.json"

# Read and parse the JSON file
const_dict = JSON.parsefile(json_file_path)

# Extract constants and their details from the dictionary
const U7_c = const_dict["U7_c"]["val"]
const U7_units = const_dict["U7_c"]["units"]
const U7_desc = const_dict["U7_c"]["desc"]

const U8_c = const_dict["U8_c"]["val"]
const U8_units = const_dict["U8_c"]["units"]
const U8_desc = const_dict["U8_c"]["desc"]

# Print the constants and their details to verify
println("U7_c = ", U7_c, " units: ", U7_units, " description: ", U7_desc)
println("U8_c = ", U8_c, " units: ", U8_units, " description: ", U8_desc)

# Define the additional function f(t)
function f(t)
    return sin(t)  # Example function, replace with your actual function
end

# Define the right-hand side function of the ODE
function rhs!(y, dy, p, t)
    f_t = f(t)  # Compute the function value at t
    dy[1] = U7_c * sin(t) - U8_c * y[1] + f_t
    dy[2] = U7_c * cos(t) - U8_c * y[2] + f_t
    dy[3] = exp(-t) - U7_c * y[3] + f_t
    dy[4] = 1/(1+t^2) - y[4]^3 + f_t
    return dy
end

# Initial conditions
y0 = [0.5, 1.0, -0.5, 0.1]

# Time span: 24 hours in seconds
t_span = (0.0, 24 * 3600.0)  # 24 hours in seconds

# Define the ODE problem
prob = ODEProblem(rhs!, y0, t_span)

# Define time evaluation points: every 10 minutes (600 seconds)
t_eval = 0:300:(24 * 3600)  # Every 600 seconds (10 minutes)

# Measure the time taken to solve the ODE using BenchmarkTools
@btime solve($prob, Tsit5(), saveat=$t_eval)

# Plot the solution
sol = solve(prob, Tsit5(), saveat=t_eval)
plot(sol, xlabel="Time (t)", ylabel="Solution (y)", title="Numerical Solution of 4-Dimensional Non-Autonomous ODE", legend=:topright)
