include("parameterLoader.jl")
include("psychrometrics.jl")
using .psychrometrics

"""
    climate_model_rhs(state, A, LAI, inputs, controls)

    Computes the right-hand side of the system of ordinary differential 
    equations (ODEs) governing the climate dynamics within a greenhouse environment.

    ### Parameters:
    - `state::Tuple{Float64, Float64, Float64, Float64}`: A tuple representing the current state of the system:
        - `T1`: Canopy temperature (°K).
        - `T2`: Air temperature (°K).
        - `C1`: CO2 concentration in the greenhouse air (ppm).
        - `V1`: Vapor pressure in the greenhouse air (kPa).

    - `inputs::Vector{Float64}`: A vector of 11 inputs that represent environmental factors affecting the greenhouse:
  - `I1`: Laef area index 
  - `I2`: Total PAR radiation outside the greenhouse (W/m²).
  - `I3`: Temperature of the heating pipes (°K).
  - `I4`: Temperature of the sky (°K).
  - `I5`: Outside air temperature (°K).
  - `I6`: Temperature of the mechanical cooling system (°K).
  - `I7`: Floor temperature (°K).
  - `I8`: Wind speed (m/s).
  - `I9`: Total radiation above the canopy (W/m²).
  - `I10`: External CO2 concentration (ppm).
  - `I11`: External vapor pressure (Pa).

    - `A::Float64`: Assimilate value representing the rate of CO2 absorption by the plants (µmol/m²/s).
    - `LAI::Float64`: Leaf area index (m^2 leaf)/(m^2 greenhouse)

    - `controls::Vector{Float64}`: A vector of 12 control variables influencing system behavior:
  - Controls for ventilation, cooling, heating, etc., such as `U1`, `U2`, ..., `U12`.

    ### Outputs:
    Returns a tuple of four values representing the rates of change of the state variables:
        - `dT1`: Rate of change of the canopy temperature (°K/s).
        - `dT2`: Rate of change of the air temperature (°K/s).
        - `dC1`: Rate of change of CO2 concentration (ppm/s).
        - `dV1`: Rate of change of vapor pressure (kPa/s).

    ### Description:    
        The function calculates the time derivatives of the greenhouse's state variables, modeling:
            1. **Radiative and heat exchanges**: Between the canopy, greenhouse air, and other surfaces (e.g., floor, sky).
            2. **Vapor and CO2 fluxes**: Through processes such as stomatal resistance, mechanical ventilation, and forced air circulation.
            3. **Climate control mechanisms**: Including heating, cooling, and ventilation adjustments based on control inputs.

        The function captures key interactions that regulate greenhouse conditions, supporting the simulation 
        of climate dynamics for energy efficiency, plant growth, and environmental control.

    ### Notes:
        - This function is part of a larger system modeling approach for greenhouse climate control.
        - It involves complex interactions of heat, moisture, and CO2 exchange, requiring accurate parameterization for effective simulation.

    ### Example call
        state = (300.0, 280.0, 0.5, 1.2)  # T1, T2, C1, V1
        inputs = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
        A = 0.3
        LAI = 2.0
        controls = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]

        dT1, dT2, dC1, dV1 = climate_model_rhs(state, A, LAI,inputs,  controls)    
        
"""
function climate_model_rhs(
    state::Tuple{Float64, Float64, Float64, Float64},   # (T1, T2, C1, V1)
    A::Float64,                                         # Assimilate
    LAI::Float64,                                       # Leaf area index
    inputs::Vector{Float64},                            # [I2, ..., I11]
    controls::Vector{Float64}                           # [U1, U2, ..., U12]
    ) :: Tuple{Float64, Float64, Float64, Float64}

    # Unpack state variables
    T1, T2, C1, V1 = state
    
    # Unpack inputs
    I2, I3, I4, I5, I6, I7, I8, I9, I10, I11 = inputs
    
    # Unpack control variables
    U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12 = controls
   
    # Global Quantities
    # Approximate VP_sat (saturation vapor pressure) of thermal screen based on average temp between air (T2) and outside temp (I5)
    T_ThmScr = 0.5 * (I5 + T2)
    VPsat_ThmScr = psychrometrics.p_sat(T_ThmScr)
    
    # Greenhouse air saturation vapor pressure (VP_sat) calculated using air temperature (T2)
    VPsat_Air = psychrometrics.p_sat(T2)

    # Outside air saturation vapor pressure (VP_sat) calculated using outside temperature (I5)
    VPsat_Out = psychrometrics.p_sat(I5)

    # VP_Can: Vapour pressure of canopy
    q2 = psychrometrics.p_sat(T1)

    # Auxiliary functions (q7, q8, q9, q10) for stomatic resistance calculations
    q7 = (1 + exp(1 / gamma5 * (I9 - delta1)))^-1
    q8 = delta4 + (delta5 - delta4) * q7
    q9 = delta6 + (delta7 - delta6) * q7
    q10 = (I9 + delta2) / (I9 + delta3)

    # VP_MechCool: Vapour pressure of mechanical cooling system
    q6 = psychrometrics.p_sat(I6)

    # rf(VP_can-VP_Air): Estomatic resistance levels large deviations of Canopy and Air vapour pressure
    q5 = 1 + q9 * ((q2 - V1)^2)

    # rf(CO2_Air): Estomatic resistance levels by large CO2 values
    q4 = 1 + q8 * ((C1 - 200)^2)

    # Estomatic resistance Canopy
    q3 = gamma4 * q10 * q4 * q5

    # VEC_CanAir: Vapour pressure exchange coefficient between Air and Canopy
    q1 = (2 * rho3 * alpha5 * LAI) / (gamma * gamma2 * (gamma3 + q3))

    # Ventilation rate equations
    n1 = nu1 * (1 - eta10 * U5)   # C_d: Discharge coefficient
    n2 = nu3 * U6                 # A_U_side: Opening of side windows
    n3 = nu2 * (1 - eta11 * U5)   # C_w: Wind pressure coefficient

    # f''_VentSide: Ventilation rate with roof vents closed
    f5 = n1 * n2 * I8 * sqrt(n3) / (2 * alpha6)
    
    # Leakage heat rate in windows         
    if I8 < 0.25
        f6 = 0.25 * nu4           # f_Leakage: Leakage rate 
    else
        f6 = nu4 * I8
    end

    # Leakage CO2 rate in windows  
    if I8 < 0.25
        f6CO2 = 0.25 * nu4CO2           # f_Leakage: Leakage rate 
    else
        f6CO2 = nu4CO2 * I8
    end

    # f''_VentRoof: Ventilation rate due to roof
    f7t = max(omega1 * nu6 * (T2 - I5) / (T2 + I5) + n3 * I8^2, 0.0)
    f7 = (U8 * nu5 * n1) / (2.0 * alpha6) * sqrt(f7t) * sign(T2 - I5)

    # Total ventilation rates due to side and roof windows, forced ventilation, and leakage
    if eta7 >= eta8
        f2 = eta6 * f5 + 0*0.5 * f6  # f_VentSide: (Total) Ventilation rate due to side windows
    else
        f2 = eta6 * (U1 * f5) + 0*0.5 * f6
    end

    f3 = U7 * phi8 / alpha6        # f_VentForced: Rate of forced ventilation
    f4 = eta6 * f7 + 0*0.5 * f6      # f_VentRoof: (Total) ventilation rate of the roof windows

    # PAR and FIR radiation calculations
    g1 = 0.49 * (1 - exp(-beta3 * LAI))               # F_PipeCan: Canopy view factor for heating pipes
    r4 = (1 - eta1) * tau1 * eta2 * I2               # R_PAR_Gh: PAR radiation above canopy
    r2 = r4 * (1 - rho1) * (1 - exp(-beta1 * LAI))    # R_PAR_SunCan_down: Direct PAR absorbed by canopy
    r3 = r4 * exp(-beta1 * LAI) * rho2 * (1 - rho1) * (1 - exp(-beta2 * LAI))
                                                     # R_PAR_SunCan_up: PAR reflected by floor and absorbed by canopy

    # FIR radiation calculations
    g2 = 1 - U1 * (1 - tau3)                         # F_CanSky: FIR exchange factor between canopy and sky
    g3 = tau1 * g2 * (1 - 0.49 * pi * gamma1 * phi1) * exp(-beta3 * LAI)
                                                     # F_FirSky: FIR exchange factor between canopy and floor

    # FIR exchanges
    SkyFIRBoltzmann = epsil3 * sigma * (T2^4 - I4^4)
    r11 = epsil4 * g3 * epsil3 * sigma * (I7^4 - I4^4)  # R_FlrSky: FIR floor-to-sky radiation
    r12 = epsil5 * (U1 * tau2) * SkyFIRBoltzmann        # R_ThScrSky: FIR fluxes between thermal screen and sky
    r13 = epsil6 * SkyFIRBoltzmann                      # R_Cov_eSky: FIR exchange between cover and sky
    r14 = epsil4 * epsil2 * (1 - 0.49 * pi * gamma1 * phi1) * sigma * (T2^4 - I7^4)
                                                      # R_CanFlr: FIR exchange between canopy and floor
    r15 = epsil2 * epsil5 * (1 - exp(-beta2 * LAI)) * g2 * sigma * (T1^4 - I4^4)
                                                      # R_CANCov_in: FIR exchange between canopy and cover

    # dT1: Canopy temperature flux
    r1 = r2 + r3                                      # R_PAR_SunCan: PAR absorbed by canopy
    r5 = (1 - eta1) * alpha2 * (eta3 * I2 + eta14 * alpha12 * U12)
                                                      # R_NIR_SunCan: NIR absorbed by canopy
    r6 = alpha3 * epsil1 * epsil2 * g1 * sigma * (I3^4 - T1^4)
                                                      # R_PipeCan: FIR exchange between canopy and pipes
    h1 = 2 * alpha4 * LAI * (T1 - T2)                  # H_CanAir: Sensible heat exchange between canopy and air
    l1 = gamma2 * q1 * (q2 - V1)                      # L_CanAir: Latent heat flux caused by transpiration
    r7 = (g1 / 0.49) * epsil2 * epsil3 * g2 * sigma * (T1^4 - I4^4)
                                                      # R_CanSky: FIR exchange between canopy and sky

    dT1 = (1.0 / (alpha1 * LAI)) * (r1 + r5 + r6 - h1 - l1 - r7 - r14)

    # dT2: Air temperature flux

    # H{PadAir} Sensible heat is exchanged between the greenhouse air and the outlet air of a cooling pad
    h2 = (U2 * phi7 / alpha6) * (rho3 * alpha5 * I5 - gamma2 * rho3 * (eta5 * (phi5 - phi6))) 

    # H{MechAir} Sensible heat is exchanged between the greenhouse air and the outlet air of a mechanical cooling system
    h3 = (U3 * lamb1 * lamb2 / alpha6) * (I6 - T2) / (T2 - I6 + (6.4e-9 * gamma2 * (V1 - q6)))  

    # H{PipeAir} Sensible heat is exchanged between the greenhouse air and the heating pipes
    h4 = n_pipes * (1.99 * pi * phi1 * gamma1 * abs((I3 - T2))^0.32) * (I3 - T2)

    # H{AirOut} Sensible heat is exchanged between the greenhouse air and the outdoor air
    # the 0.5*f6 comes form separation of CO2 and heat leakage
    h7 = rho3 * alpha5 * (f2 + f3 + 0.5*f6) * (T2 - I5)

    # H{AirFlr} Sensible heat is exchanged between the greenhouse air and the floor
    h11 = 2 * nu7 * (T2 - I7) / (phi2 + nu8)

    # Sensible heat due to lamps
    h12 = (eta15 + eta16) * alpha12 * U12

    # R{GlobSunAir} Global radiation absorbed by construction elements and released to the air
    r8 = I2 * (eta1 * (tau1 * eta2 + (alpha2 + alpha7) * eta3) + (alpha8 * eta2 + alpha9 * eta3))

    # R{AirSky} Total radiation to the sky
    r10 = r11 + r12 + r13

    # L{AirFog} Latent heat needed to evaporate water droplets from the fogging system
    l2 = gamma2 * U9 * phi9 / alpha6

    # Calculate the temperature derivative of T2
    dT2 = (1 / (phi2 * rho3 * alpha5)) * (h1 + h2 + h4 - h7 - h11 + h12 + r8 - r10 - l2)
   
    # dV1: Vapor preasure flux

    # MV{CanAir} Vapour exchanged between air and canopy
    p1 = q1 * (q2 - V1)

    # MV{PadAir} Vapour exchanged between air and the outlet air of the pad
    p2 = rho3 * (U2 * phi7 / alpha6) * (eta5 * (phi5 - phi6) + phi6)

    # MV{FogAir} Vapour exchanged between air and fogging system
    p3 = U9 * phi9 / alpha6

    # MV{BlowAir} Vapour exchanged between air and the direct air heater
    p4 = eta12 * U4 * lamb4 / alpha6

    # MV{AirOut} Vapour exchanged between the greenhouse air and outdoor air
    p5 = (psi1 / omega2) * (V1 / T2 - I11 / I5) * (f2 + f3 + f4 + f6)

    # MV{AirOutPad} Vapour exchanged between the air and the outdoor air due to the air exchange from the pad and fan system
    p6 = (U2 * phi7 / alpha6) * (psi1 / omega2) * (V1 / T2)

    # MV{AirMech} Vapour exchanged between air and the mechanical cooling system
    p7 = V1 < q6 ? 0 : 6.4e-9 * (V1 - q6)

    # MV{AirThScr} Vapour exchanged between air and the thermal screen (not used in this model)
    p8 = V1 < VPsat_ThmScr ? 0 : 6.4e-9 * 1.7 * (abs(T2 - VPsat_ThmScr))^0.33 * (V1 - VPsat_ThmScr)

    # MV{AirTop} Vapour exchanged between air and the top
    p9 = V1 < VPsat_Air ? 0 : 6.4e-9 * 1.2 * 1000 * ((abs(T1 - T2))^0.33) * (V1 - VPsat_Air)

    # Calculate the vapour derivative of V1
    kappa3 = (psi1 * phi2) / (omega2 * T2)
    dV1 = (1 / kappa3) * (p1 + p2 + p3 + p4 - p5 - p6 - p7 -  p8 + p9)

    # dC1: CO2 mass flux 

    # MC{BlowAir} CO2 flux from the heat blower to the greenhouse air
    o1 = eta13 * U4 * lamb4 / alpha6

    # MC{ExtAir} CO2 exchanged between the greenhouse air and external CO2 source
    o2 = U10 * psi2 / alpha6

    # MC{PadAir} CO2 exchanged between the greenhouse air and the pad and fan system
    o3 = (U2 * phi7 / alpha6) * (C1 - I10)

    # MC{AirCan} CO2 flux between the greenhouse air and the canopy
    o4 = A

    # MC{AirOut} CO2 flux between the greenhouse air and outside air
    o5 = (C1 - I10) * (f2 + f3 + f4 + f6CO2)

    # Calculate the CO2 derivative of C1
    dC1 = (1 / phi2) * (o1 + o2 + o3 - o4 - o5)
   
    return   dT1, dT2, dC1, dV1
end


"""
state = (300.0, 280.0, 0.5, 1.2)  # T1, T2, C1, V1
inputs = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
A = 0.3
LAI = 2.0
controls = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]

dT1, dT2, dC1, dV1 = climate_model_rhs(state, A, LAI,inputs,  controls)   
"""