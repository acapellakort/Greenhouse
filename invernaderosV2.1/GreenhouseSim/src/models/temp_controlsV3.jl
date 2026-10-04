# Temporary definitions of control functions for system testing

# Control function for the thermal screen (currently inactive)
function U1f(t)
    return 0.0  # Screen is fully open
end

# Control for the fan-pad system(currently inactive)
function U2f(t)
    return 0.0  # No effect
end

# Control for the Mechanical cooling system (currently inactive)
function U3f(t)
    return 0.0  # Fan-pad system is off
end

# Control for the air heater (currently inactive)
function U4f(t)
    return 0.0  # Heater is off
end

# Control for the shading screen (currently inactive)
function U5f(t)
    return 0.0  # Screen is fully open
end

# Side windows opening (currently inactive)
function U6f(t)
    return 0.0  # No effect
end

# Control for forced ventilation (partially active)
function U7f(t)
    return 0.0  # Ventilation at 50% capacity
end

# Roof windows control (partially active)
function U8f(t)
    return 0.0  # Placeholder, 50% effect
end

# Fog System control (partially active)
function U9f(t)
    return 0.5  # Placeholder, 50% effect
end

# Control for the CO2 source (currently inactive)
function U10f(t)
    return 0.0  # CO2 source is off
end

# Control for pipe temperature (active)
function U11f(t)
    return 23 + 273.15  # Pipe temperature set to 23°C (296.15 K)
end

# Light control  (active)
function U12f(t)
    return 1.0  #
end
