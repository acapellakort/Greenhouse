using DifferentialEquations
using Plots

# Define parameters
T_set = 35.0        # Setpoint temperature (°C)
T_ambient = 20.0     # Ambient temperature (°C)
tau = 10.0           # Thermal time constant of the chamber
Kp = 2.0             # Proportional gain
Ki = 0.1             # Integral gain
decay_factor = 0.01  # Decay factor for the integral term (small value to slowly decay old errors)
t_final = 100.0      # Simulation duration

MaxPot = 100        # Maximum power of the heatiung system

# Initial conditions
T0 = T_ambient       # Initial temperature (°C)
integral_error0 = 0.0  # Initial integral of the error

# Define the PI control function with leaky integral
function pi_control(T, integral_error, Kp, Ki, decay_factor, dt)
    # Compute the error
    error = T_set - T
    
    # Update the integral of error with decay (leaky integral)
    integral_error = integral_error * (1 - decay_factor * dt) + error * dt
    
    # PI control law
    P = Kp * error + Ki * integral_error
    
    # Clamp the heater power to between 0 and 100%
    P = clamp(P, 0, MaxPot)
    
    return P, integral_error
end

# ODE system: [T(t), integral_error]
function pi_control_ode!(du, u, p, t)
    T, integral_error = u  # Unpack the current state variables (Temperature, Integral Error)

    # Calculate a constant small dt (Euler step approximation) decay rate
    dt = 0.01
    
    # Calculate control signal (heater power) using the PI controller with leaky integral
    P, integral_error = pi_control(T, integral_error, Kp, Ki, decay_factor, dt)

    # ODE for the temperature
    du[1] = (P - (T - T_ambient)) / tau 
    # ODE for the integral of the error (for the integral part of the PI control)
    du[2] = T_set - T  # This is the error that will accumulate for the next step

    u[2] = integral_error  # Update integral_error for the next iteration
end

# Initial state: temperature and integral of error
u0 = [T0, integral_error0]

# Time span
tspan = (0.0, t_final)

# Define the ODE problem (parameters p are not used, so pass nothing)
prob = ODEProblem(pi_control_ode!, u0, tspan, nothing)

# Solve the ODE problem
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

# Extract temperature and power values for plotting
T_vals = sol[1, :]        # Temperature values
integral_error_vals = sol[2, :]  # Integral of error (for reference)
P_vals = pi_control.(T_vals, integral_error_vals, Kp, Ki, decay_factor, 0.01)[1]  # Compute power using the PI control law

# Plot the results
plot(sol.t, T_vals)
