
# Function for K_C
function K_C(T)
    pow = 0.1*(T - 298.15) 
    return K_C25 * Q10_KC^pow
end

# Function for K_O
function K_O(T)
    pow = 0.1*(T - 298.15) 
    return K_O25 * Q10_KO^pow
end


########################## Model for stomatal conductance 

# Function for r_s (stomatal resistance)
function r_s(T, I, C, V, LAI)
    # T: Temperature (K)
    # I: Incident light (W/m^2)
    # C: CO2 concentration (mg/m^3)
    # V: Vapor pressure (Pa)
    # LAI: Leaf area index (dimensionless)
    
    C_i = 0.509e-6 * C  # Internal CO2 concentration (mol/m^3)
    I2 = I_2(I)  # Light intensity transformed by a specific function (W/m^2)
    
    # Saturation ratio Sr
    Sr = 1 / (1 + exp(S * (I2 - Rs)))  # Sr is a sigmoid function based on I2 and Rs (dimensionless)
    
    # CO2 resistance factor
    C_ev3 = C_ev3n * (1 - Sr) + C_ev3d * Sr  # Interpolates between two constants C_ev3n and C_ev3d (dimensionless)
    C_ev4 = C_ev4n * (1 - Sr) + C_ev4d * Sr  # Interpolates between two constants C_ev4n and C_ev4d (Pa^-1)
    
    f_R = (I2 / LAI + C_ev1) / (I2 / LAI + C_ev2)  # Radiation function (dimensionless)
    
    # CO2 correction factor
    f_C = 1 + C_ev3 * (C_i - 200e-6)^2  # (dimensionless)
    
    VP_sat = Pws(T)  # Saturation vapor pressure at temperature T (Pa)
    
    rh = rhf(T, VP_sat)  # Relative humidity based on T and VP_sat (%)
    
    PD = VPDf(V, rh)  # Vapor pressure deficit (Pa)
    
    f_V = 1 + C_ev4 * PD^2  # Vapor pressure correction factor (dimensionless)

    # Returns stomatal resistance in units of s/m
    return r_m * f_R * f_C * f_V
end

# Simple alternative model for stomatal conductance
function g_eff()
    # Based on Bush (2023)
    # gsc: Stomatal conductance (mol/m²/s)
    gsc  = 0.25  # Typical value
    
    # gm: Mesophyll conductance (mol/m²/s)
    gm   = 0.3  # Typical value

    # Returns effective conductance (mol/m²/s)
    return gm * gsc / (gm + gsc)
end

################## Quadratic functions with resistance to calculate assimilates

function AcropFast(T, C, Iw, LAI)
    # T: Temperature (K)
    # C: CO2 concentration (mg/m^3)
    # Iw: Incident light (W/m^2)
    # LAI: Leaf area index (dimensionless)
    
    # Convert incident light from W/m^2 to mol/m²/s
    I = 4.6e-6 * Iw
    
    # Convert CO2 concentration to mol/m³
    C_i = 0.509e-6 * C
    
    # Gamma_st_t: Compensation point for CO2 (mol/m^3)
    Gamma_st_t = (0.5 * O_a) / (Sco25 * exp(((T - 298.15) / T) * (E_Soc) / (298.15 * Rgas)))
    
    # Electron transport rate (mol/m²/s)
    J_t = ((alpha * I + J_max) - sqrt((alpha * I + J_max)^2 - 4 * theta * alpha * I * J_max)) / (2 * theta)
    
    # Effective conductance (mol/m²/s)
    g = g_eff()
    
    # Maximum carboxylation rate (mol/m²/s)
    V_cmax_t = V_cmax25 * Q10_Vcmax^(0.1 * (T - 298.15)) / (1 + exp(0.128 * (T - 315.15)))
    
    # Michaelis-Menten constants for CO2 and O2 (mol/m³)
    K_C_t = K_C(T)
    K_O_t = K_O(T)

    # Calculate assimilation rates
    # For RuBisCO-limited photosynthesis
    pc = -(V_cmax_t + g * (C_i + K_C_t * (1 + O_a / K_O_t)) - Rd_day)
    qc = g * (V_cmax_t * max(C_i - Gamma_st_t, 0) - (C_i + K_C_t * (1 + O_a / K_O_t)) * Rd_day)
    disc_c = pc^2 - 4 * qc
    Ar_c = disc_c >= 0 ? 0.5 * (-pc - sqrt(disc_c)) : -Rd_day

    # For light-limited photosynthesis
    pj = -0.25 * J_t + Rd_day - g * (C_i + 2.0 * Gamma_st_t)
    qj = 0.25 * g * max(C_i - Gamma_st_t, 0) * J_t - g * (C_i + 2 * Gamma_st_t) * Rd_day
    disc_j = pj^2 - 4 * qj
    Ar_j = disc_j >= 0 ? 0.5 * (-pj - sqrt(disc_j)) : -Rd_day

    # For product-limited photosynthesis
    Ar_p = 1.5 * V_cmax_t - Rd_day

    # Calculate final assimilation rate (mg/m²/s)
    # Convert mol/m²/s to mg/m²/s using molar mass of CO2 (44.010 mg/mol)
    return 44010 * LAI * minimum([Ar_j, Ar_c, Ar_p])
end
