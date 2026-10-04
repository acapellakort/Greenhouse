# Function to calculate V_cmax
function V_cmax(T)
    #pow1 = 0.1*(T - 298.15) 
    #pow2 = 0.128 * (T - 315.15) 
    return (V_cmax25) * Q10_Vcmax^(0.1*(T - 298.15)) / (1 + exp(0.128 * (T - 315.15)))
end

function Gamma_st(T)
    # Gamma_st based on Xin & Struick 2009
    # Units in paper are of mol_CO2 * mol_Air^-1 
    pow = ((T - 298.15)/T) * (E_Soc) / (298.15 * Rgas)
    return (0.5 * O_a)/ (Sco25 * exp(pow))
end

# Tau function
function tau(T)
    pow = (T - 298.15) / 10
    return tau_25 * Q10_tau^pow
end

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

# Function for I_2
function I_2(I)
    return alpha*I 
end

# Function for J
function J(I)
    # [I] = mol_photons m^-2 s^-1   
 return ((alpha*I + J_max) - sqrt((alpha*I + J_max)^2 - 4 * theta * alpha*I * J_max)) / (2 * theta)
end



# Function W_C according to Bush 2023
function W_C(T, C)
    # C enters in units of mg * m^-3 we trasnform it to 1e6ppm ie unitless 
    C_c = 0.509e-6 * C
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    Wc = V_cmax(T)*C_c/(C_c + K_C(T)*(1 + O_a/K_O(T)) )
    return Wc
end

# Function W_J according to Bush 2023
function W_J(T, C, I)
    # C enters in units of mg * m^-3
    C_c = 0.506e-6 * C
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    WJ =  J(I) * C_c/ (4 * C_c + 8 * Gamma_st(T))
    return WJ
end


# Function W_p according to Bush 2023
function W_p(T, C)
    # C enters in units of mg * m^-3
    C_c = 0.509e-6 * C
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    Tpmax = V_cmax(T)/2.0
    Wp =  3*Tpmax * maximum([C_c/ (C_c - Gamma_st(T)), 1])
    return Wp
end



# according to Bush 2023
function A_B(T,C,I)
    # C enters in units of mg * m^-3
    C_c = 0.509e-6 * C
    factor = minimum([1-Gamma_st(T)/C_c,1])
    A = minimum([W_C(T,C),W_J(T,C,I),W_p(T,C)]) * factor - Rd_day
    return A
end


# Function for A_acum
function A_acum(T)
    return V_cmax(T)/2.0
end

# Function for R_d
function R_d(V_cmax)
    return 0.015 * V_cmax
end


# Function for Acrop
function Acrop(A, L)
    return (30 / 10^6) * L * A
end

function Pws(T1)
    # saturation water vapor pressure measure in Pa
    q2 = ifelse.(T1 .> 273.15,
        611.21 .* exp.((18.678 .- ((T1 .- 273.15) ./ 234.5)) .* ((T1 .- 273.15) ./ (257.14 .+ T1 .- 273.15))),
        611.15 .* exp.((23.036 .- ((T1 .- 273.15) ./ 333.7)) .* ((T1 .- 273.15) ./ (279.82 .+ T1 .- 273.15)))
    )
    return q2
end

# Relatiuve Humidity
function rhf(T1, V1)
    # Usando saturation water vapor pressure measure in Pa y VP en pascales
    return 100 .* V1 ./ Pws(T1)
end

# Vapour pressure deficit
function VPDf(VP_sat, rh) 
    return VP_sat.*(1 .- rh./100)
end

########################## modelo para conductancia estomatica 

# Function for r_s (stomatal resistance)
function r_s(T, I, C, V, LAI )
    # C enters in units of mg * m^-3
    C_i = 0.509e-6 * C
    I2 = I_2(I)
    # Function for Sr (saturation ratio)
    Sr = 1 / (1 + exp(S * (I2 - Rs)))
    # Function for C_ev3 (CO2 resistance factor calculation)
    C_ev3 = C_ev3n * (1 - Sr) + C_ev3d * Sr #[-]
    C_ev4 = C_ev4n*(1 -Sr) + C_ev4d*Sr #[Pa^-1]
    f_R = (I2 / LAI + C_ev1) / (I2 / LAI + C_ev2) # [-]
    f_C = 1 + C_ev3 * (C_i - 200e-6)^2  # [-]
    VP_sat = Pws(T) #[Pa]
    rh = rhf(T, VP_sat) #[P%]
    PD = VPDf(V, rh)
    f_V = 1 + C_ev4 * PD^2  #  [-]

     #[r_m] s /m 
    return r_m * f_R * f_C * f_V
end

function gsf( rs )
    # rs es la resistencia estomatica calculada antes
    # toma la conductancia del agua en unidade de  s * m**-1
    # la conductancia del CO2 es (1/1.6) veces la del agua
    # regresa la resistencia en unidades de 
    #  mu_mol_CH2O * m**-1 * s**-1 * ppm_CO2
    # Los factores de conversion son:
    # 1/rs *(0.553/0.044)
    return (0.553/0.044)*(1.6 * rs)^-1
end

# Resistencia estomática para el agua en unidades de s * m**-1
# R_agua = r_s( r_m=self.V('r_m'), f_R=f_R1, f_C=f_C1, f_V=f_V1) 


# alternative simple model for stomatal conductnace
function g_eff()
    # Bush 2023
    # Stomatal conductance
    gsc  = 0.25 
    # mesophyll conductance 
    gm   = 0.3 

    return gm*gsc/(gm+gsc)
end

################## Funciones cuadratcias con resistencia para calcualr Assimilates

###############################################################


function Ar_J(T,C,I)
    C_i = 0.509e-6 * C
     # The resuls is in mol * m^-2 * s^-1 that makes sense
    Gamma_st_t = Gamma_st(T) 
    g = g_eff()
    p = -0.25*J(I) + Rd_day - g * (C_i + 2.0 * Gamma_st_t) 
    q = 0.25*g*(max(C_i - Gamma_st_t, 0)) * J(I) - g*(C_i + 2 * Gamma_st_t)*Rd_day
    if p^2 - 4*q >= 0
        return 0.5*(-p - sqrt(p^2 - 4*q))  
    else
        return -Rd_day
    end
end

function Ar_C(T,C)
    C_i =0.509e-6 * C
     # The resuls is in mol * m^-2 * s^-1 that makes sense
    Gamma_st_t = Gamma_st(T)
    g = g_eff()
    Vcmax_t = V_cmax(T)
    K_C_t = K_C25 * Q10_KC^(0.1*(T - 298.15))
    K_O_t = K_O25 * Q10_KO^(0.1*(T - 298.15))
    p = -(Vcmax_t + g*(C_i + K_C_t*(1 + O_a/K_O_t))- Rd_day)
    q = g *(Vcmax_t*max(C_i - Gamma_st_t, 0) - (C_i + K_C_t*(1 + O_a/K_O_t))* Rd_day ) 
    if p^2 - 4*q >= 0
        return 0.5*(-p - sqrt(p^2 - 4*q))  
    else
        return -Rd_day
    end
end


# Function W_p according to Bush 2023
function Ar_P(T)
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    Wp =  1.5*V_cmax(T) - Rd_day
    return Wp
end

function Aresult(T,C,I)
    # C enters in units of mg * m^-3
    A = minimum([Ar_J(T,C,I),Ar_C(T,C),Ar_P(T)]) 
    return A
end

function Acrop(T, C, Iw, LAI)
    # T en K ok
    # C en mg * m^-3 ok
    # Iw esta en Watts por metro^2 y lo necesitamos en mol m⁻² s⁻¹
    I = 4.6e-6 * Iw
    Ar = Aresult(T,C,I)
    # Para el modele necesitamos A en "mg * m^-2 * s^-1"
    # donde m^2 son de invernadero. 
    # A result nos da el resulatado en "mol m⁻² s⁻¹", pero 
    # en m^2 de hoja. Por lo tango, tenemos que trasnformar
    #
    # 1 mol de CO2 son 44,010 mg de CO2
    #
    return 44010 * Ar*LAI
end


function AcropFast(T, C, Iw, LAI)
    # T en K ok
    # C en mg * m^-3 ok
    # Iw esta en Watts por metro^2 y lo necesitamos en mol m⁻² s⁻¹
    I = 4.6e-6 * Iw
    C_i = 0.509e-6 * C
     # The resuls is in mol * m^-2 * s^-1 that makes sense
     # Global quantities
    Gamma_st_t = (0.5 * O_a)/ (Sco25 * exp(((T - 298.15)/T) * (E_Soc) / (298.15 * Rgas)))
    J_t = ((alpha*I + J_max) - sqrt((alpha*I + J_max)^2 - 4 * theta * alpha*I * J_max)) / (2 * theta)
    g = g_eff()
    V_cmax_t = (V_cmax25) * Q10_Vcmax^(0.1*(T - 298.15)) / (1 + exp(0.128 * (T - 315.15)))
    K_C_t = K_C(T)
    K_O_t = K_O(T)

    # Calculo Ar_C
    pc = -(V_cmax_t + g*(C_i + K_C_t*(1 + O_a/K_O_t))- Rd_day)
    qc = g *(V_cmax_t*max(C_i - Gamma_st_t, 0) - (C_i + K_C_t*(1 + O_a/K_O_t))* Rd_day ) 
    disc_c = pc^2 - 4*qc
    if disc_c >= 0
        Ar_c = 0.5*(-pc - sqrt(disc_c))  
    else
        Ar_c = -Rd_day
    end

    # Calculo Ar_J
    pj = -0.25*J_t + Rd_day - g * (C_i + 2.0 * Gamma_st_t) 
    qj = 0.25*g*(max(C_i - Gamma_st_t, 0)) * J_t - g*(C_i + 2 * Gamma_st_t)*Rd_day
    disc_j = pj^2 - 4*qj 
    if disc_j >= 0
        Ar_j = 0.5*(-pj - sqrt(disc_j))  
    else
        Ar_j = -Rd_day
    end

    # Calculo Ar_P
    Ar_p =  1.5*V_cmax_t - Rd_day

    # Ar = minimum([Ar_J,Ar_c,Ar_P]) 
    # Para el modele necesitamos A en "mg * m^-2 * s^-1"
    # donde m^2 son de invernadero. 
    # A result nos da el resulatado en "mol m⁻² s⁻¹", pero 
    # en m^2 de hoja. Por lo tango, tenemos que trasnformar
    #
    # 1 mol de CO2 son 44,010 mg de CO2
    #
    return 44010 * LAI * minimum([Ar_j, Ar_c, Ar_p])
end