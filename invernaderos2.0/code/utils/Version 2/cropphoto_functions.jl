# Function to calculate V_cmax
function V_cmax(T)
    pow1 = (T - 298.15u"K") / (10u"K")
    pow2 = 0.128 * (T - 315.15u"K") / (1u"K")
    return (V_cmax25) * Q10_Vcmax^(pow1) / (1 + exp(pow2))
end

function Gamma_st(T)
    # Gamma_st based on Xin & Struick 2009
    # Units in paper are of mol_CO2 * mol_Air^-1 
    pow = ((T - 298.15u"K")/T) * (E_Soc) / (298.15u"K" * Rgas)
    return (0.5 * O_a)/ (Sco25 * exp(pow))
end

# Tau function
function tau(T)
    pow = (T - 298.15u"K") / 10u"K"
    return tau_25 * Q10_tau^pow
end

# Function for K_C
function K_C(T)
    pow = (T - 298.15u"K") / (10u"K")
    return K_C25 * Q10_KC^pow
end

# Function for K_O
function K_O(T)
    pow = (T - 298.15u"K") / (10u"K")
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
    C_c = 0.509e-6u"mg^-1 * m^3"* C
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    Wc = V_cmax(T)*C_c/(C_c + K_C(T)*(1 + O_a/K_O(T)) )
    return Wc
end

# Function W_J according to Bush 2023
function W_J(T, C, I)
    # C enters in units of mg * m^-3
    C_c = 0.506e-6u"mg^-1 * m^3"* C
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    WJ =  J(I) * C_c/ (4 * C_c + 8 * Gamma_st(T))
    return WJ
end


# Function W_p according to Bush 2023
function W_p(T, C)
    # C enters in units of mg * m^-3
    C_c = 0.509e-6u"mg^-1 * m^3"* C
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    Tpmax = V_cmax(T)/2.0
    Wp =  3*Tpmax * maximum([C_c/ (C_c - Gamma_st(T)), 1])
    return Wp
end



# according to Bush 2023
function A_B(T,C,I)
    # C enters in units of mg * m^-3
    C_c = 0.509e-6u"mg^-1 * m^3"* C
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

# Function for A
function A(A_R, A_f, A_acum)
    return min(A_R, A_f, A_acum)
end

# Function for Acrop
function Acrop(A, L)
    return (30 / 10^6) * L * A
end

function Pws(T1)
    # saturation water vapor pressure measure in Pa
    q2 = ifelse.(T1 .> 273.15u"K",
        611.21u"Pa" .* exp.((18.678 .- ((T1 .- 273.15u"K") ./ 234.5u"K")) .* ((T1 .- 273.15u"K") ./ (257.14u"K" .+ T1 .- 273.15u"K"))),
        611.15u"Pa" .* exp.((23.036 .- ((T1 .- 273.15u"K") ./ 333.7u"K")) .* ((T1 .- 273.15u"K") ./ (279.82u"K" .+ T1 .- 273.15u"K")))
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
    C_i = 0.556e-6u"mg^-1 * m^3"* C
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
    gsc  = 0.25u"mol * m^-2*s^-1"
    # mesophyll conductance 
    gm   = 0.3u"mol * m^-2*s^-1"

    return gm*gsc/(gm+gsc)
end

################## Funciones cuadratcias con resistencia para calcualr Assimilates

# Function for Ar_f
function Ar_f(T,C,VP,I,LAI)
    # C enters in units of mg * m^-3
    C_i = 0.556e6u"mg^-1 * m^3"* C
     # The resuls is in mol * m^-2 * s^-1 that makes sense
    Gamma_st_t = Gamma_st(T)
    rs = r_s(T, I, C, VP, LAI)
    g = g_eff()
    p = 0.25*J(I) + Rd_day + g * (C_i - 2.0 * Gamma_st_t) 
    q = 0.25*g*(max(C_i - Gamma_st_t, 0)) * J(I) + g*(C_i + 2 * Gamma_st_t)*Rd_day

    return (-p - sqrt(p^2 - 4*q)) / 2
    
end

###############################################################


function Ar_J(T,C,I)
    C_i = 0.509e-6u"mg^-1 * m^3"* C
     # The resuls is in mol * m^-2 * s^-1 that makes sense
    Gamma_st_t = Gamma_st(T) 
    g = g_eff()
    p = -0.25*J(I) + Rd_day - g * (C_i + 2.0 * Gamma_st_t) 
    q = 0.25*g*(max(C_i - Gamma_st_t, 0)) * J(I) - g*(C_i + 2 * Gamma_st_t)*Rd_day
    if p^2 - 4*q >= 0
        return (-p - sqrt(p^2 - 4*q)) / 2
    else
        return -Rd_day
    end
end

function Ar_C(T,C)
    C_i = 0.509e-6u"mg^-1 * m^3"* C
     # The resuls is in mol * m^-2 * s^-1 that makes sense
    Gamma_st_t = Gamma_st(T)
    g = g_eff()
    Vcmax_t = V_cmax(T)
    K_C_t = K_C(T)
    K_O_t = K_O(T)
    p = -(Vcmax_t + g*(C_i + K_C_t*(1 + O_a/K_O_t))- Rd_day)
    q = g *(Vcmax_t*max(C_i - Gamma_st_t, 0) - (C_i + K_C_t*(1 + O_a/K_O_t))* Rd_day ) 
    if p^2 - 4*q >= 0
        return (-p - sqrt(p^2 - 4*q)) / 2
    else
        return -Rd_day
    end
end


# Function W_p according to Bush 2023
function Ar_P(T)
    # The resuls is in mol * m^-2 * s^-1 that makes sense
    Wp =  3*V_cmax(T)/2.0 -Rd_day
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
    I = 4.6e-6u"mol * m^-2 * s^-1 * m^2 * W^-1" * Iw
    Ar = Aresult(T,C,I)
    # Para el modele necesitamos A en "mg * m^-2 * s^-1"
    # donde m^2 son de invernadero. 
    # A result nos da el resulatado en "mol m⁻² s⁻¹", pero 
    # en m^2 de hoja. Por lo tango, tenemos que trasnformar
    #
    # 1 mol de CO2 son 44,010 mg de CO2
    #
    return 44010u"mg* mol^-1"*Ar*LAI
end