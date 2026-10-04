
# Models for Stomatal Conductance #

"""

    r_s(T::Float64, I::Float64, C::Float64, V::Float64, LAI::Float64) -> Float64

#  Compute stomatal resistance depending on envorinmental variables 

# Arguments:
- T::Float64: Temperature in Kelvin
- I::Float64: Incident light intensity (W/m²)
- C::Float64: CO₂ concentration (mg/m³)
- V::Float64: Vapor pressure (Pa)
- LAI::Float64: Leaf Area Index (dimensionless)

# Returns:
- Float64: Stomatal resistance (s/m)
"""
function r_s(T::Float64, I::Float64, C::Float64, V::Float64, LAI::Float64)::Float64
    C_i = 0.509e-6 * C  # Internal CO₂ concentration (mol/m³)
    I2 = I_2(I)         # Transformed light intensity (W/m²)
    
    # Saturation ratio, calculated as a sigmoid function
    Sr = 1 / (1 + exp(S * (I2 - Rs)))
    
    # Interpolating CO₂ resistance factors
    C_ev3 = C_ev3n * (1 - Sr) + C_ev3d * Sr  # Dimensionless
    C_ev4 = C_ev4n * (1 - Sr) + C_ev4d * Sr  # Pa^-1
    
    # Radiation correction factor
    f_R = (I2 / LAI + C_ev1) / (I2 / LAI + C_ev2)
    
    # CO₂ correction factor
    f_C = 1 + C_ev3 * (C_i - 200e-6)^2  # Dimensionless
    
    # Calculate vapor pressure and related corrections
    VP_sat = Pws(T)    # Saturation vapor pressure (Pa)
    rh = rhf(T, VP_sat)  # Relative humidity (%)
    PD = VPDf(V, rh)     # Vapor pressure deficit (Pa)
    f_V = 1 + C_ev4 * PD^2  # Dimensionless

    # Return stomatal resistance (s/m)
    return r_m * f_R * f_C * f_V
end

"""
    g_eff()::Float64

Compute effective conductance combining stomatal and mesophyll conductance
in a simplified model based on  Bush (2023).

# Arguments: None

# Returns:
- Float64: Effective conductance (mol/m²/s)
"""
function g_eff()::Float64
    gsc = 0.25  # Stomatal conductance (mol/m²/s)
    gm = 0.3    # Mesophyll conductance (mol/m²/s)
    
    # Effective conductance using harmonic mean
    return gm * gsc / (gm + gsc)
end

# Models for CO₂ and O₂ Michaelis-Menten constants #
"""

    K_C(T::Float64)->Float64

Compute Michaelis-Menten constant for CO₂ at temperature `T`

# Arguments:
- T::Float64: Temperature in Kelvin

# Returns:
- Float64: Michaelis-Menten constant for CO₂ (mol/m³)
"""
function K_C(T::Float64)::Float64
    pow = 0.1 * (T - 298.15)
    return K_C25 * Q10_KC^pow
end

"""

    K_O(T::Float64)->Float64

Compute Michaelis-Menten constant for O₂ at temperature `T`

# Arguments:
- T::Float64: Temperature in Kelvin

# Returns:
- Float64: Michaelis-Menten constant for CO₂ (mol/m³)
"""
function K_O(T::Float64)::Float64
    pow = 0.1 * (T - 298.15)
    return K_O25 * Q10_KO^pow
end

# Photosynthesis model #

"""
    compute_assimilates(T::Float64, C::Float64, Iw::Float64, LAI::Float64) -> Float64

Calculate the rate of assimilates in the chloroplast using 
the FvCB (Farquhar-von Caemmerer-Berry) photosynthesis model. 
This function also integrates an effective conductance model 
that accounts for both stomatal and mesophyll conductance.

# Arguments
- `T::Float64`: Canopy temperature in Kelvin.
- `C::Float64`: CO₂ concentration in mg/m³.
- `Iw::Float64`: Incident light intensity at the crop in W/m².
- `LAI::Float64`: Leaf Area Index (dimensionless).

# Returns
A `Float64` value representing the assimilate production rate.

# Example call

    T = 300.0 
    C = 400.0 
    Iw = 1000.0
    LAI = 2.0  
    assimilates = compute_assimilates(T, C, Iw, LAI)
"""
function compute_assimilates(T::Float64, C::Float64, Iw::Float64, LAI::Float64)::Float64

    # Convert incident light (W/m²) to moles of photons (mol/m²/s)
    I = 4.6e-6 * Iw
    
    # Convert CO₂ concentration (mg/m³) to mol/m³
    C_i = 0.509e-6 * C
    
    # Gamma_st_t: CO₂ compensation point (mol/m³)
    Gamma_st_t = (0.5 * O_a) / (Sco25 * exp(((T - 298.15) / T) * (E_Soc) / (298.15 * Rgas)))
    
    # Electron transport rate (mol/m²/s) using the light response curve
    J_t = ((alpha * I + J_max) - sqrt((alpha * I + J_max)^2 - 4 * theta * alpha * I * J_max)) / (2 * theta)
    
    # Effective conductance (mol/m²/s) combining stomatal and mesophyll conductance
    g = g_eff()
    
    # Maximum carboxylation rate (mol/m²/s), adjusted for temperature
    V_cmax_t = V_cmax25 * Q10_Vcmax^(0.1 * (T - 298.15)) / (1 + exp(0.128 * (T - 315.15)))
    
    # Michaelis-Menten constants for CO₂ and O₂ (mol/m³), dependent on temperature
    K_C_t = K_C(T)
    K_O_t = K_O(T)

    # Calculate RuBisCO-limited photosynthesis rate (Ar_c)
    # pc and qc are coefficients for solving the quadratic assimilation rate equation
    pc = -(V_cmax_t + g * (C_i + K_C_t * (1 + O_a / K_O_t)) - Rd_day)
    qc = g * (V_cmax_t * max(C_i - Gamma_st_t, 0) - (C_i + K_C_t * (1 + O_a / K_O_t)) * Rd_day)
    disc_c = pc^2 - 4 * qc
    Ar_c = disc_c >= 0 ? 0.5 * (-pc - sqrt(disc_c)) : -Rd_day

    # Calculate light-limited photosynthesis rate (Ar_j)
    # pj and qj are coefficients for solving the quadratic assimilation rate equation
    pj = -0.25 * J_t + Rd_day - g * (C_i + 2.0 * Gamma_st_t)
    qj = 0.25 * g * max(C_i - Gamma_st_t, 0) * J_t - g * (C_i + 2 * Gamma_st_t) * Rd_day
    disc_j = pj^2 - 4 * qj
    Ar_j = disc_j >= 0 ? 0.5 * (-pj - sqrt(disc_j)) : -Rd_day

    # Calculate product-limited photosynthesis rate (Ar_p)
    Ar_p = 1.5 * V_cmax_t - Rd_day

    # Convert the final assimilation rate from mol/m²/s to mg/m²/s
    # Using the molar mass of CO₂ (44.010 mg/mol)
    return 44010 * LAI * minimum([Ar_j, Ar_c, Ar_p])
end


