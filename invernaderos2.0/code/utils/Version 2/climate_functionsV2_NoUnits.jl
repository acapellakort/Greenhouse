
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

######################
# Functions for T1 rhs

function Rr1(I1, I2)
    r4 = (1 - eta1) * tau1 * eta2 .* I2
    r2 = r4 .* (1 - rho1) .* (1 .- exp(-beta1 .* I1))
    r3 = r4 .* exp(-beta1 .* I1) * rho2 * (1 - rho1) .* (1 .- exp(-beta2 .* I1))
    return r2 + r3
end

function Rr5(I2, U12)
    return (1 - eta1) .* alpha2 * (eta3 .* I2 .+ eta14 * alpha12 .* U12)
end

function Rr6(T1, I1, I3)
    g1 = 0.49 .* (1 .- exp(-beta3.* I1))
    return alpha3 * epsil1 * epsil2 * g1 .* sigma * (I3.^4 - T1.^4)
end

function Rh1(T1, T2, I1)
    return 2 * alpha4 * I1 .* (T1 .- T2)
end

function Rl1(T1, C1, V1, I1, I9)
    q10 = (I9 .+ delta2) ./ (I9 .+ delta3)
    q7 = (1 .+ exp.(gamma5 .* (I9 .- delta1))).^-1
    q9 = delta6 .* (1 .- q7) + delta7 .* q7
    q8 = delta4 .* (1 .- q7) + delta5 .* q7

    q2 = ifelse.(T1 .> 273.15,
    611.21 .* exp.((18.678 .- ((T1 .- 273.15) ./ 234.5)) .* ((T1 .- 273.15) ./ (257.14 .+ T1 .- 273.15))),
    611.15 .* exp.((23.036 .- ((T1 .- 273.15) ./ 333.7)) .* ((T1 .- 273.15) ./ (279.82 .+ T1 .- 273.15)))
    )

    q5 = 1 .+ q9 .* ((q2 .- V1).^2)
    q4 = 1 .+ q8 .* ((eta4 .* C1 .- 200).^2)
    q3 = gamma4 .* q10 .* q4 .* q5
    q1 = (2 .* rho3 .* alpha5 .* I1) ./ (gamma * gamma2 .* (gamma3 .+ q3))
    p1 = q1 .* (q2 .- V1)
    
    return gamma2 .* p1
end

function Rr7(T1, I1, I4, U1)
    a1 = 1 .- exp(-beta3 .* I1)
    b1 = 1 .- U1 * (1 - tau3)
    g2 = tau2 .* b1
    return a1 .* epsil2 * epsil3 .* g2 * sigma .* (T1.^4 .- I4.^4)
end

function Rr14(T2, I7) #- R_{CanFlr}
    return epsil4 .* epsil2 .*  (1 .- 0.49 .* pi .* gamma1 .* phi1) .* sigma .* (T2.^4 .-I7.^4) 
end

function Rr15(T1, I1, I4, U1)  #R_{CANCov,in}
    g2 = 1 .- U1 .* (1 .- tau3)  
    return  epsil2 .* epsil5 .* (1 .-exp.(-beta2.*I1)).* g2.* sigma .* (T1.^4 .-  I4.^4) 
end


# Functions for T2 rhs

function Rh2(I5, U2) 
    return (U2 * phi7 / alpha6) * (rho3 * alpha5 * I5 - gamma2 * rho3 * (eta5 * (phi5 - phi6)))
end 

function Rh3(T2, V1, I6, U3)

    if I6 > 273.15
       q6 = 611.21 * exp((18.678 - ((I6-273.15) / 234.5)) * ((I6 - 273.15)/ (257.14 + I6-273.15)))
    else
       q6 = 611.15 * exp((23.036 - ((I6-273.15) / 333.7)) * ((I6-273.15)/ (279.82 + I6-273.15)))
    end

    num = (U3 * lamb1 * lamb2 / alpha6) * ( I6 - T2) 
    # The units in  denominator's constant are NOT in the Thesis, but needed for consistency
    den = T2 - I6 + (6.4e-9 * gamma2 * (V1 - q6))

    return num/den
end  

function Rh4(T2, I3)
    HEC = n_pipes*(1.99 * pi * phi1 * gamma1 * abs((I3 - T2))^0.32 )
    return  HEC *   (I3 - T2)
end

function Rh5(T2, I7)
    return sigma3 * (I7 - T2)
end 

function Rh7(T2, I5, I8, U1, U5, U6, U7)

    f3 = U7 * phi8 / alpha6
    f5 = nu1 * (1 - eta10 * U5) * (nu3 * U6) * I8 * sqrt(nu2*(1 - eta11*U5)) / (2 * alpha6)

    if I8 < 0.25 
        f6 = 0.25 * nu4 
    else 
        f6 = nu4 * I8
    end

    if eta7 >= eta8 
        f2 = eta6 * f5 + 0.5 * f6 
    else 
        f2 = eta6 * (U1 * f5) + 0.5 * f6
    end 

    return  rho3 * alpha5 * (f2 + f3) * (T2 - I5)
end

function Rh11(T2, I7)
    return 2 * nu7 * (T2 - I7)/(phi2 + nu8)
end

function Rh12(U12)
    return (eta15 + eta16 ) * alpha12 * U12
end

function Rr10(T2, I1, I4, U1)
    g3 = tau1 * (1 - U1 * (1 - tau3)) * (1 - 0.49 * pi * gamma1 * phi1) * exp(-beta3 * I1)
    r11 = epsil4 * epsil3 * g3 * sigma * (T2^4 -  I4^4)
    r12 = epsil5 * epsil3 * (U1 * tau2) * sigma * ( T2^4 -  I4^4 )
    r13 = epsil6 * epsil3 * sigma * ( T2 ^4 -  I4^4)
    return r11 + r12 + r13
end

function Rr8(I2) 
    return eta1 * I2 * (tau1 * eta2 + (alpha2 + alpha7)*eta3) + (alpha8 * eta2 + alpha9 * eta3) * I2

end

function Rl2(U9)
    return gamma2 * U9 * phi9 / alpha6
end

# Functions for V1 rhs
function Rp1(T1, C1, V1, I1, I9)

    if T1 > 273.15
        q2 = 611.21 * exp((18.678 - ((T1-273.15) / 234.5)) * ((T1 - 273.15)/ (257.14 + T1-273.15)))
    else
        q2 = 611.15* exp((23.036 - ((T1-273.15) / 333.7)) * ((T1-273.15)/ (279.82 + T1-273.15)))
    end

    q8  = delta4  + (delta5 - delta4) * (1 + exp((1/gamma5) * (I9 - delta1)))^-1
    q9  = delta6  + (delta7 - delta6)* (1 + exp((1/gamma5) * (I9 - delta1)))^-1
    q10 = (I9 + delta2) / (I9 + delta3)
    q3  = gamma4 * q10 * (1 + q8 * ((eta4 * C1 - 200)^2)) * (1 + q9 * ((q2 - V1)^2)) 
    q1  = (2 * rho3 * alpha5 * I1) / (gamma * gamma2 * (gamma3 + q3))
    return q1 * (q2 - V1)
end 

function Rp2(U2)
    return rho3 * (U2 * phi7 / alpha6) * (eta5 * (phi5 - phi6) + phi6)
end 

function Rp3(U9)
    return U9 * phi9 / alpha6
end    

function Rp4(U4)
    return eta12 * U4 * lamb4 / alpha6
end

function Rp5(T2, V1, I5, I8, I11, U1, U5, U6, U7, U8 )

    n1 = nu1*(1 - eta10 * U5)
    n2 = nu3 * U6 
    n3 = nu2 * (1 - eta11 * U5)

    f5 = n1 * n2 * I8 * sqrt(n3) / (2 * alpha6)
    if I8 < 0.25 
        f6 = 0.25 * nu4 
    else 
        f6 = nu4 * I8
    end

    f7t = omega1 * nu6 * (T2 - I5) / (T2 + I5)+ n3 * I8^2
    f7 = (U8 * nu5 * n1) / (2.0 * alpha6) * sqrt(max( f7t, 0.0))
  
    if eta7 >= eta8
        f2 =  eta6 * f5 + 0.5 * f6 
    else
        f2 = eta6 * (U1 * f5) + 0.5 * f6
    end
    f3 = U7 * phi8 / alpha6
    f4 = eta6 * f7 + 0.5 * f6 
    
   
    return (psi1 / omega2) * (V1 / T2 - I11 / I5) * (f2 + f3 + f4)

end 

function Rp6(T2, V1, U2)
    return (U2 * phi7 / alpha6) * (psi1 / omega2) * (V1 / T2)
end

function Rp7(T2, V1, I6, U3)
    h3 = Rh3(T2, V1, I6,  U3)

    if I6 > 273.15
        q6 = 611.21 * exp((18.678 - ((I6-273.15) / 234.5)) * ((I6 - 273.15)/ (257.14 + I6-273.15)))
     else
        q6 = 611.15 * exp((23.036 - ((I6-273.15) / 333.7)) * ((I6-273.15)/ (279.82 + I6-273.15)))
    end
 
 
    if V1 < q6 
        return 0 
    else  
        return 6.4e-9 * h3 * (V1 - q6)
    end
end

# Functions for C1 rhs

function Ro1(U4) 
    return eta13 * U4*lamb4/alpha6
end

function Ro2(U10) 
    return U10 * psi2 / alpha6
end

function Ro3(C1, I10, U2) 
    return (U2 * phi7 / alpha6) * (C1 - I10)
end 

function Ro4(A) 
    return A
end

function Ro5(T2, C1, I5, I8, I10, U1, U5, U6, U7, U8) 
    n1 = nu1*(1 - eta10 * U5)
    n2 = nu3 * U6 
    n3 = nu2 * (1 - eta11 * U5)

    f5 = n1 * n2 * I8 * sqrt(n3) / (2 * alpha6)
    if I8 < 0.25 
        f6 = 0.25 * nu4 
    else 
        f6 = nu4 * I8
    end

    f7 = omega1 * nu6 * (T2 - I5) / (T2 + I5)+ n3 * I8^2
    f7 = (U8 * nu5 * n1) / (2.0 * alpha6) * sqrt(max( f7, 0.0))
  
    if eta7 >= eta8
        f2 =  eta6 * f5 + 0.5 * f6 
    else
        f2 = eta6 * (U1 * f5) + 0.5 * f6
    end
    f3 = U7 * phi8 / alpha6
    f4 = eta6 * f7 + 0.5 * f6 
 
    return  (C1 - I10) * (f2 + f3 + f4)
end 

  
# rhs functions ----------------------------------------------

function rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12)
    # Flux functions as they appear in the Thesis
    r_1 = Rr1(I1, I2)
    r_5 = Rr5(I2, U12)
    r_6 = Rr6(T1, I1, I3)
    h_1 = Rh1(T1, T2, I1)
    l_1 = Rl1(T1, C1, V1, I1, I9)
    r_7 = Rr7(T1, I1, I4, U1)

    return (1.0/(alpha1 * I1))*(r_1 + 0.2*r_5 + r_6 - h_1 - l_1 - r_7)
end

function rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)
    # Flux functions as they appear in the Thesis
    h_1  = Rh1(T1, T2, I1) 
    h_2  = Rh2(I5, U2) 
    #h_3  = Rh3(T2, V1, I6, U3) +h3 pero es cero y las unidades estan mal. Las formulas no estan claras
    h_4  = Rh4(T2, I3)

    h_7  = Rh7(T2, I5, I8, U1, U5, U6, U7)
    h_11 = Rh11(T2, I7)
    h_12 = Rh12(U12)
    r_8  = Rr8(I2)
    r_10 = Rr10(T2, I1, I4, U1)
    l_2 = Rl2(U9)

    #return (1/(phi2 * rho3 * alpha5))*(h_1 + h_2 + h_4 + h_5 - h_7 - h_11 + h_12 + r_8 - r_10  - l_2 )
    return (1/(phi2 * rho3 * alpha5))*(h_1 + h_2 + h_4 - h_7 - h_11 + h_12 + r_8 - r_10  - l_2 )

end

function rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7, U8, U9)
    # Flux functions as they appear in the Thesis
    p_1  = Rp1(T1, C1, V1, I1, I9) 
    p_2  = Rp2(U2) 
    p_3  = Rp3(U9)
    p_4  = Rp4(U4)
    p_5  = Rp5(T2, V1, I5, I8, I11, U1, U5, U6, U7, U8)
    p_6  = Rp6(T2, V1, U2)
    p_7  = Rp7(T2, V1, I6, U3)
    kappa3 = (psi1 * phi2) / (omega2 * T2)

    return (1/kappa3) * (p_1 + p_2 + p_3 + p_4 - p_5 - p_6 - p_7)    
end

function rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10)
    # Flux functions as they appear in the Thesis
    o_1  = Ro1(U4)  
    o_2  = Ro2(U10) 
    o_3  = Ro3(C1, I10, U2)
    o_4  = Ro4(A) 
    o_5  = Ro5(T2, C1, I5, I8, I10, U1, U5, U6, U7, U8)
    return (1/phi2)*(o_1 + o_2 + o_3 - o_4 - o_5 )  
end


# main rhs ----------------------------------------------

# Define the right-hand side function of the ODE
function rhs!(dy, y, p, t)
    # In this function y = [T1, T2, C1, V1])
        #Tout = I5
        #WindSpeed = I8
        #CO2out = I10 
        #Iglobal = I2   
        #VPout = I11
        # I1 = LAI  
        # I2 = Iglobal  
        # I3 = Tpipef    
        # I4 = Tskif
        # I5 = Tout
        # I6 = TmechCool  
        # I7 = Tsoil 
        # I8 = WindSpeed m/s
        # I9 = Idocel
        #I10 = CO2out  
        #I11 = VPout

    
    # Move evaluation of weather to current date
    t_floatR = t  + p
    #weather data
    I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11 = 
     LAIf(t_floatR), Iglobalf(t_floatR),  Tpipef(t_floatR), 
     Tskyf(t_floatR), Toutf(t_floatR), TmechCoolf(t_floatR), 
     Tsoilf(t_floatR), WindSpeedf(t_floatR), Idocelf(t_floatR),  
     CO2outf(t_floatR),  VPoutf(t_floatR)

    # Control variables this is a temporary implementation no real time
    U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12 =  
    U1f(t),  U2f(t), U3f(t), U4f(t), U5f(t), U6f(t), 
    U7f(t), U8f(t), U9f(t), U10f(t), U12f(t)
    
    #Assimilates function photosynthesis   
    A = Af(t)    
    
    # unpack y 
    T1, T2, C1, V1 = y[1], y[2], y[3], y[4]

     # compute RHS
    dy[1] = rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12) #u"K * s^-1" 
    dy[2] = rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)#u"K * s^-1" 
    dy[3] = rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10) # u"mg * m^-3 * s^-1"
    dy[4] = rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7,U8, U9)#u"Pa* s^-1"

    #return [dy1, dy2, dy3, dy4]
end


    



############### Fast version for the functionsfunction rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7, U8, U9)
    # Flux functions as they appear in the Thesis
   
function rhs_fast( T1, T2, C1, V1, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A, U1,  U2, U3,  U4, U5, U6, U7,U8, U9,U10, U12)

    # Global quantities

    # Approximation of VP_sat of the thermal screen as the one for the average temp 
    # between air and oiutside temps
    T_ThmScr = 0.5*(I5+T2)
    VPsat_ThmScr = ifelse.(T_ThmScr .> 273.15,
    611.21 .* exp.((18.678 .- ((T_ThmScr .- 273.15) ./ 234.5)) .* ((T_ThmScr .- 273.15) ./ (257.14 .+ T_ThmScr .- 273.15))),
    611.15 .* exp.((23.036 .- ((T_ThmScr .- 273.15) ./ 333.7)) .* ((T_ThmScr .- 273.15) ./ (279.82 .+ T_ThmScr .- 273.15)))
    )
    # Green house air saturation vapor pressure VPsat_air
    VPsat_air = ifelse.(T2 .> 273.15,
    611.21 .* exp.((18.678 .- ((T2 .- 273.15) ./ 234.5)) .* ((T2 .- 273.15) ./ (257.14 .+ T2 .- 273.15))),
    611.15 .* exp.((23.036 .- ((T2 .- 273.15) ./ 333.7)) .* ((T2 .- 273.15) ./ (279.82 .+ T2 .- 273.15)))
    )

    # Ouside air saturation preasure
    VPsat_out = ifelse.(I5 .> 273.15,
    611.21 .* exp.((18.678 .- ((I5 .- 273.15) ./ 234.5)) .* ((I5 .- 273.15) ./ (257.14 .+ I5 .- 273.15))),
    611.15 .* exp.((23.036 .- ((I5 .- 273.15) ./ 333.7)) .* ((I5 .- 273.15) ./ (279.82 .+ I5 .- 273.15)))
    )

    # q2 = VP_Can Vapour pressure canopy
    if T1 > 273.15
        q2 = 611.21 * exp((18.678 - ((T1-273.15) / 234.5)) * ((T1 - 273.15)/ (257.14 + T1-273.15)))
    else
        q2 = 611.15* exp((23.036 - ((T1-273.15) / 333.7)) * ((T1-273.15)/ (279.82 + T1-273.15)))
    end
    
    q7 = (1 + exp(1/gamma5 * (I9 - delta1)))^-1     # q7 = Auxiliary function
    q8  = delta4  + (delta5 - delta4) * q7          # q8 = Auxiliary function
    q9  = delta6  + (delta7 - delta6)* q7           # q9 = Auxiliary function
    q10 = (I9 + delta2) / (I9 + delta3)             # q10= Auxiliary function

    # VP_MechCool Vapour pressure  Mechanical Cooling System    
    if I6 > 273.15                                  # VP_MachCool Vapour pressure  Mechanical Cooling System
        q6 = 611.21 * exp((18.678 - ((I6-273.15) / 234.5)) * ((I6 - 273.15)/ (257.14 + I6-273.15)))
     else
        q6 = 611.15 * exp((23.036 - ((I6-273.15) / 333.7)) * ((I6-273.15)/ (279.82 + I6-273.15)))
    end

    q5 = 1 + q9 * ((q2 - V1)^2)                     # rf(VP_can-VP_Air) Estomatic resistance levels large deviations of Canopy and Air vapour pressure
    q4 = 1 + q8 * ((eta4 * C1 - 200)^2)             # rf(CO_{2Air-ppm}) Estomatic resistance levels by large CO2 values
    q3  = gamma4 * q10 * q4 * q5                    # Estomatic resistance Canopy
    q1  = (2 * rho3 * alpha5 * I1) / (gamma * gamma2 * (gamma3 + q3)) 
                                                    # VEC_{CanAir} VP interchange coefficient bewteen Air and Canopy 
    #--------------------------------------
    n1 = nu1*(1 - eta10 * U5)                       # C_d discharge coefficient
    n2 = nu3 * U6                                   # A^U_side Opening of side windows 
    n3 = nu2 * (1 - eta11 * U5)                     # C_w Global coefficient of preasure due to WindSpeed

    f5 = n1 * n2 * I8 * sqrt(n3) / (2 * alpha6)     #f''_{VentSide} Ventilation rate with roof vents closed
    if I8 < 0.25 
        f6 = 0.25 * nu4                             #f_{Leakage} Leakage rate 
    else 
        f6 = nu4 * I8
    end

    f7t = max(omega1 * nu6 * (T2 - I5) / (T2 + I5)+ n3 * I8^2, 0.0)
    f7 = (U8 * nu5 * n1) / (2.0 * alpha6) * sqrt(f7t )*sign(T2-I5)  # f''_{VentRoof} Ventilation rate due to roof
  
    if eta7 >= eta8
        f2 =  eta6 * f5 + 0.5 * f6                  # f_{VentSide} Total ventilation rate due to side windows
    else
        f2 = eta6 * (U1 * f5) + 0.5 * f6
    end
    f3 = U7 * phi8 / alpha6                         # f_{VentForced} Rate of forced venntilation
    f4 = eta6 * f7 + 0.5 * f6                       # f_{VentRoof} Tptal ventilation rate of the roof windows
    
    # Radiacion 
    # PAR
    g1 =  0.49*(1 - exp(-beta3 * I1))                # F_{PipeCan} Factor de vista        
    r4 = (1 - eta1) * tau1 * eta2 * I2               # R_{Par_Gh} PAR radiation avobe canopy
                                                     # Eq 27
    r2 = r4 * (1 - rho1) * (1 - exp(-beta1 * I1))    # R_{PAR_SunCan}down down PAR transmitted by the greenhouse cover that is directly absorbed by the canopy
                                                     # Eq 26
    r3 = r4 * exp(-beta1 * I1) * rho2 * (1 - rho1) * (1 - exp(-beta2 * I1))
                                                     # R_{PAR_SunCan}up PAR that is first reflected by the greenhouse floor and then is absorbed by the canopy
                                                     # Eq 28

    # FIR radiation 
    g2 = 1 - U1 * (1 - tau3)                         # F_{CanSky} Factor FIR is exchanged between the canopy  and the sky 
    g3 = tau1 * g2 * (1 - 0.49 * pi * gamma1 * phi1) * exp(-beta3 * I1) 
                                                     # F_{FirSky} Factor FIR is exchanged between the canopy  and the floor
    SkyFIRBoltzmann  = epsil3 * sigma * (T2^4 -  I4^4) 
    r11 = epsil4 * g3 * epsil3 * sigma * (I7^4 -  I4^4) # R_{FlrSky} FIR floor radiation to sky 
    r12 = epsil5 * (U1 * tau2) * SkyFIRBoltzmann     # R_{ThScrSky} FIR fluxes between the thermal screen and  sky  
    r13 = epsil6 * SkyFIRBoltzmann                   # R_{Cov,eSky} FIR exchange between the cover and the sky.
    r_14 =  epsil4 * epsil2 *  (1 - 0.49 * pi * gamma1 * phi1) * sigma * (T2^4-I7^4) 
                                                      #- R_{CanFlr}
    r_15 =  epsil2 * epsil5 * (1-exp(-beta2*I1))*g2* sigma * (T1^4 -  I4^4) 
                                                      #- R_{CANCov,in}
    #Inicio: Calculo de T1 #############################

    r_1 = r2 + r3                                   # R_{PAR_SunCan} The PAR absorbed by the canopy  
                                                    # eq (25) Vanthoor
    r_5 = (1 - eta1) * alpha2 * (eta3 * I2 + eta14 * alpha12 * U12) 
                                                    # R_{NIR_SunCan} The NIR absorbed by the canopy
                                                    # Eq 32
    r_6 = alpha3 * epsil1 * epsil2 * g1 * sigma * (I3^4 - T1^4)
                                                    # R_{PipeCan} FIR is exchanged between the canopy and heating pipes 
    h_1  = 2 *alpha4 * I1 * (T1 - T2)              # H_{CanAir} is the sensible heat exchange between canopy and greenhouse air 
                                                    #Eq 39
    l_1 = gamma2 * q1 * (q2 - V1)                   # L_{CanAir} is the latent heat flux caused by transpiration.
    r_7 = (g1/0.49) * epsil2 * epsil3 * g2 * sigma * (T1^4 - I4^4)
                                                    # R_{CanSky} FIR is exchanged between the canopy and sky
    

    #T1 = (1.0/(alpha1 * I1))*(r_1 + 0.2*r_5 + r_6 - h_1 - l_1 - r_7)
    dT1= (1.0/(alpha1 * I1))*(r_1 + r_5 + r_6 - h_1 - l_1 - r_7-  r_14- 0*r_15)

    # R_{PAR_SunCan} r_1 + R_{NIR_SunCan} r_5 + R_{PipeCan} r_6
    #        - H_{CanAir} h_1 - L_{CANAir} l_1 - R_{CANCov,in} - R_{CanFlr} r_14 - R_{CanSky} r_7
    #Fin: Calculo de T1 ###################

    #Inicio: Calculo de T2 #############################
    
    h_2  = (U2 * phi7 / alpha6) * (rho3 * alpha5 * I5 - gamma2 * rho3 * (eta5 * (phi5 - phi6)))
                                                    # H_{PadAir} Sensible heat is exchanged between the greenhouse air and the outlet air of a cooling pad  <-----------
    h_3  = (U3 * lamb1 * lamb2 / alpha6) * ( I6 - T2) /( T2 - I6 + (6.4e-9 * gamma2 * (V1 - q6)))
                                                    # H_{MechAir} Sensible heat is exchanged between the greenhouse air and the outlet air of a Mechanical cooling system 
    h_4 =  n_pipes*(1.99 * pi * phi1 * gamma1 * abs((I3 - T2))^0.32 ) *   (I3 - T2)
                                                    # H_{PipeAir} Sensible heat is exchanged between the greenhouse air and the outlet air of a Mechanical cooling system the heating pipes 
    #h_5 = sigma3 * (I7 - T2) sigma3 not defined    # H_{PasAir} Sensible heat is exchanged between the greenhouse air and the passive energy buffer
    h_7 = rho3 * alpha5 * (f2 + f3) * (T2 - I5)     # H_{AirOut} Sensible heat is exchanged between the greenhouse air and the outdoor air
    h_11 = 2 * nu7 * (T2 - I7)/(phi2 + nu8)         # H_{AirFlr} Sensible heat is exchanged between the greenhouse air and the floor
    h_12 = (eta15 + eta16 ) * alpha12 * U12         # Sensible heat due to lamps
    r_8  = I2 *(eta1 * (tau1 * eta2 + (alpha2 + alpha7)*eta3) + (alpha8 * eta2 + alpha9 * eta3) )    
                                                    # R_{Glob_SunAir} global radiation which is absorbed by the construction elements and which is released to the air
    r_10 = r11 + r12 + r13                          # R_{AirSky}   total radiation to the sky not in Vanthoor as such
    l_2 = gamma2 * U9 * phi9 / alpha6               # L_{AirFog}  latent heat needed to evaporate the water droplets added by a fogging system.

    #T3 =  (1/(phi2 * rho3 * alpha5))*(h_1 + h_2 + h_4 + h_5 - h_7 - h_11 + h_12 + r_8 - r_10  - l_2 )
    dT2 = (1/(phi2 * rho3 * alpha5))*(h_1 + h_2 + h_4 - h_7 - h_11 + h_12 + r_8 - r_10  - l_2 )
    #Fin: Calculo de T2 ###################

    #Inicio: Calculo de V1 ###################
    p_1 =  q1 * (q2 - V1)                           # MV_{CanAir} Vapour is exchanged between the air and the canopy
    p_2 = rho3 * (U2 * phi7 / alpha6) * (eta5 * (phi5 - phi6) + phi6)
                                                    # MV_{PadAir} Vapour is exchanged between the air and the outlet air of the pad 
    p_3 = U9 * phi9 / alpha6                        # MV_{FogAir} Vapour is exchanged between the air and  the fogging system
    p_4 = eta12 * U4 * lamb4 / alpha6               # MV_{BlowAir} Vapour is exchanged between the air and the direct air heater 
    p_5 = (psi1 / omega2) * (V1 / T2 - I11 / I5) * (f2 + f3 + f4) 
                                                    # MV_{AirOut} Vapour is exchanged between the (greenhouse) air and  outdoor air 
    p_6 = (U2 * phi7 / alpha6) * (psi1 / omega2) * (V1 / T2)
                                                    # MV_{AirOut_Pad} Vapour is exchanged between the air and the outdoor air due to the air exchange caused by the pad and fan system
    if V1 < q6 
        p_7 =  0                                    # MV_{AirMech} Vapour is exchanged between the air and and the mechanical cooling system 
    else  
        p_7 = 6.4e-9 * (V1 - q6)
    end
    # We add the condensation on the termal screen assuming that the 
    # termperature of the termal screen is the temperature of the outside.
    # We also assume that that the thermal screen is completely closed  
    #   - MV_{AirThScr} 
    if V1 < VPsat_ThmScr  # al final no lo usamos
        p_8 =  0                                    # MV_{AirThScr} Vapour is exchanged between the air and and the thermal screen always 
    else  
        p_8 = 6.4e-9 * 1.7*(abs(T2-VPsat_ThmScr ))^0.33*(V1 -  VPsat_ThmScr )
    end

    #- MV_{AirTop}
    #   - MV_{AirThScr} 
    #   - MV_{AirCoc,in}  No se como hacerle!!!!
    if V1 < VPsat_air
        p_9 =  0                                    # MV_{AirThScr} Vapour is exchanged between the air and and the thermal screen always 
    else  
        p_9 = 6.4e-9 * 1.2 * 1000 * ((abs(T1-T2))^0.33 ) *(V1 -  VPsat_air)
    end
    

    kappa3 = (psi1 * phi2) / (omega2 * T2)

    dV1 = (1/kappa3) * (p_1 + p_2 + p_3 + p_4 - p_5 - p_6 - p_7 - 0*p_8 - p_9)  
    # MV_{CanAir} p_1 + MV_{PadAir} p_2 +  MV_{FogAir} p_3 + MV_{BlowAir} p_4 
    #           -  MV_{AirOut} p_5 - MV_{AirOut_Pad} p_6 - MV_{AirMech}* p_7

    # Vanthoor 
    # MV_{CanAir}* p_1 + MV_{PadAir}* p_2 +  MV_{FogAir}* p_3 + MV_{BlowAir}* p_4 
    #          - MV_{AirThScr} - MV_{AirTop} 
    # -  MV_{AirOut}* p_5 - MV_{AirOut_Pad}* p_6 - MV_{AirMech} *p_7
    #
    #Fin: Calculo de V1 ###################
    
    
    #Inicio: Calculo de C1 ###################
    o_1  = eta13 * U4*lamb4/alpha6                  # MC_{BlowAir} CO2 flux from the heat blower to the greenhouse air 
    o_2  = U10 * psi2 / alpha6                      # MC_{ExtAir} Carbon dioxide is exchanged between the greenhouse air and external CO2 source
    o_3  = (U2 * phi7 / alpha6) * (C1 - I10)        # MC_{PadAir} Carbon dioxide is exchanged between the greenhouse air  and the pad and fan system 
    o_4  = A                                        # MC_{AirCan} CO2 flux between the greenhouse air and the canopy 
    o_5 = (C1 - I10) * (f2 + f3 + f4)               # MC_{AirOut} CO2 flux between the greenhouse air and outside air
    dC1 = (1/phi2)*(o_1 + o_2 + o_3 - o_4 - o_5 )  
    #Fin: Calculo de C1 ####################### 

        return   dT1, dT2, dC1, dV1
    end
    
#precompile(rhs_fast, (Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,))    