# Functions for T1 rhs

function Rr1(I1, I2)
    r2 = (1 - eta1) * tau1 * eta2 * I2 * (1 - rho1) * (1 - exp(-beta1 * I1))
    r3 = (1 - eta1) * tau1 * eta2 * I2 * exp(-beta1 * I1) * rho2 * (1 - rho1) * (1 - exp(-beta2 * I1))
    return r2 + r3
end

function Rr5(I2, U12)
    return (1 - eta1) * alpha2 * (eta3 * I2 + eta14 * alpha12 * U12)
end

function Rr6(T1, I1, I3)
    g1 = 0.49 * (1 - exp(-beta3 * I1))
    return alpha3 * epsil1 * epsil2 * g1 * sigma * (I3^4 - T1^4)
end

function Rh1(T1, T2, I1)
    return 2 * alpha4 * I1 * (T1 - T2)
end

function Rl1(T1, C1, V1, I1, I9)
    q10 = (I9 + delta2) / (I9 + delta3)
    q7 = (1 + exp(1/gamma5 * (I9 - delta1)))^-1
    q9 = delta6 * (1 - q7) + delta7 * q7
    q8 = delta4 * (1 - q7) + delta5 * q7

    if T1 > 0u"°C"
        q2 = 611.21u"Pa" * exp((18.678 - ((uconvert(u"K", T1)-273.15u"K") / 234.5u"K")) * ((uconvert(u"K", T1) -273.15u"K")/ (257.14u"K" + (T1-273.15u"K"))))
    else
        q2 = 611.15u"Pa" * exp((23.036 - ((uconvert(u"K", T1)-273.15u"K") / 333.7u"K")) * ((uconvert(u"K", T1) -273.15u"K")/ (279.82u"K" + (T1-273.15u"K"))))
    end

    q5 = 1 + q9 * ((q2 - V1)^2)
    q4 = 1 + q8 * ((eta4 * C1 - 200u"ppm")^2)
    q3 = gamma4 * q10 * q4 * q5
    q1 = (2 * rho3 * alpha5 * I1) / (gamma * gamma2 * (gamma3 + q3))
    p1 = q1 * (q2 - V1)
    
    return gamma2 * p1
end

function Rr7(T1, I1, I4, U1)
    a1 = 1 - exp(-beta3 * I1)
    b1 = 1 - U1 * (1 - tau3)
    g2 = tau2 * b1
    return a1 * epsil2 * epsil3 * g2 * sigma * (T1^4 - I4^4)
end

# Functions for T2 rhs

function Rh2(I5, U2) 
    return (U2 * phi7 / alpha6) * (rho3 * alpha5 * uconvert(u"K", I5) - gamma2 * rho3 * (eta5 * (phi5 - phi6)))
end 

function Rh3(T2, V1, I6, U3)

    if I6 > 0u"°C" 
       q6 = 611.21u"Pa"*exp((18.678 - ((uconvert(u"K", I6)-273.15u"K") / 234.5u"K")) * ((uconvert(u"K", I6) -273.15u"K")/ (257.14u"K" + (I6-273.15u"K"))))
    else
       q6 = 611.15u"Pa"*exp((23.036 - ((uconvert(u"K", I6)-273.15u"K") / 333.7u"K")) * ((uconvert(u"K", I6) -273.15u"K")/ (279.82u"K" + (I6-273.15u"K"))))
    end

    num = (U3 * lamb1 * lamb2 / alpha6) * ( I6 - T2) 
    # The units in  denominator's constant are NOT in the Thesis, but needed for consistency
    den = T2 - I6 + (6.4e-9u"K * kg * J^-1 * Pa^-1" * gamma2 * (V1 - q6))

    return num/den
end  

function Rh4(T2, I3)
    HEC = n_pipes*(1.99 * pi * phi1 * gamma1 * abs(ustrip((I3 - T2)))^0.32 )u"W * m^-2 * K^-1"
    return  HEC *  uconvert(u"K",(I3 - T2))
end

function Rh5(T2, I7)
    return sigma3 * (I7 - T2)
end 

function Rh7(T2, I5, I8, U1, U5, U6, U7)

    f3 = U7 * phi8 / alpha6
    f5 = nu1 * (1 - eta10 * U5) * (nu3 * U6) * I8 * sqrt(nu2*(1 - eta11*U5)) / (2 * alpha6)

    if I8 < 0.25u"m * s^-1"
        f6 = 0.25u"m * s^-1" * nu4 
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
    r11 = epsil4 * epsil3 * g3 * sigma * ((uconvert(u"K", T2))^4 - (uconvert(u"K", I4))^4)
    r12 = epsil5 * epsil3 * (U1 * tau2) * sigma * ((uconvert(u"K", T2))^4 - (uconvert(u"K", I4))^4)
    r13 = epsil6 * epsil3 * sigma * ((uconvert(u"K", T2))^4 - (uconvert(u"K", I4))^4)
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

    if T1 > 0u"°C"
        q2 = 611.21u"Pa" * exp((18.678 - ((uconvert(u"K", T1)-273.15u"K") / 234.5u"K")) * ((uconvert(u"K", T1) -273.15u"K")/ (257.14u"K" + (T1-273.15u"K"))))
    else
        q2 = 611.15u"Pa" * exp((23.036 - ((uconvert(u"K", T1)-273.15u"K") / 333.7u"K")) * ((uconvert(u"K", T1) -273.15u"K")/ (279.82u"K" + (T1-273.15u"K"))))
    end


    q8  = delta4  + (delta5 - delta4) * (1 + exp((1/gamma5) * (I9 - delta1)))^-1
    q9  = delta6  + (delta7 - delta6)* (1 + exp((1/gamma5) * (I9 - delta1)))^-1
    q10 = (I9 + delta2) / (I9 + delta3)
    q3  = gamma4 * q10 * (1 + q8 * ((eta4 * C1 - 200u"ppm")^2)) * (1 + q9 * ((q2 - V1)^2)) 
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
    if I8 < 0.25u"m * s^-1"
        f6 = 0.25u"m * s^-1" * nu4 
    else 
        f6 = nu4 * I8
    end

    f7 = omega1 * nu6 * (uconvert(u"K",T2) - uconvert(u"K",I5)) / (uconvert(u"K",T2) + uconvert(u"K",I5))+ n3 * I8^2
    f7 = (U8 * nu5 * n1) / (2.0 * alpha6) * sqrt(max( f7, 0.0u"m^2 * s^-2"))
  
    if eta7 >= eta8
        f2 =  eta6 * f5 + 0.5 * f6 
    else
        f2 = eta6 * (U1 * f5) + 0.5 * f6
    end
    f3 = U7 * phi8 / alpha6
    f4 = eta6 * f7 + 0.5 * f6 
    
    return (psi1 / omega2) * ((V1 / uconvert(u"K",T2)) - (I11 / uconvert(u"K",I5))) * (f2 + f3 + f4)

end 

function Rp6(T2, V1, U2)
    return (U2 * phi7 / alpha6) * (psi1 / omega2) * (V1 / uconvert(u"K",T2))
end

function Rp7(T2, V1, I6, U3)
    h3 = Rh3(T2, V1, I6,  U3)

    if I6 > 0u"°C" 
        q6 = 611.21u"Pa"*exp((18.678 - ((uconvert(u"K", I6)-273.15u"K") / 234.5u"K")) * ((uconvert(u"K", T1) -273.15u"K")/ (257.14u"K" + (I6-273.15u"K"))))
     else
        q6 = 611.15u"Pa"*exp((23.036 - ((uconvert(u"K", I6)-273.15u"K") / 333.7u"K")) * ((uconvert(u"K", I6) -273.15u"K")/ (279.82u"K" + (I6-273.15u"K"))))
     end
 
 
    if V1 < q6 
        return 0u"kg * m^-2 * s^-1"
    else  
        return (6.4e-9u"kg * Pa^-1 * W^-1 * s^-1") * h3 * (V1 - q6)
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
    if I8 < 0.25u"m * s^-1"
        f6 = 0.25u"m * s^-1" * nu4 
    else 
        f6 = nu4 * I8
    end

    f7 = omega1 * nu6 * (uconvert(u"K",T2) - uconvert(u"K",I5)) / (uconvert(u"K",T2) + uconvert(u"K",I5))+ n3 * I8^2
    f7 = (U8 * nu5 * n1) / (2.0 * alpha6) * sqrt(max( f7, 0.0u"m^2 * s^-2"))
  
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
    kappa3 = (psi1 * phi2) / (omega2 * uconvert(u"K",T2))

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
        # I3 = Tpipef   +++++++?
        # I4 = Tskif
        # I5 = Tout
        # I6 = TmechCool   =  0u"°C"
        # I7 = Tsoil =  17u"°C"
        # I8 = WindSpeed km/h
        # I9 = Idocel
        #I10 = CO2out  
        #I11 = VPout

    
    # Convert time to Float64
    t_float = t
    t_floatR = t  + p
    #weather data
    I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11 = 
     LAIf(t_floatR), Iglobalf(t_floatR),  Tpipef(t_floatR), 
     Tskyf(t_floatR), uconvert(u"K",  Toutf(t_floatR)), 
     TmechCoolf(t_floatR),  Tsoilf(t_floatR), 
     uconvert(u"m/s",  WindSpeedf(t_floatR)), 
     Idocelf(t_floatR),  CO2outf(t_floatR),  VPoutf(t_floatR)

    # Control variables
    U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12 =  
    U1f(t_float),  U2f(t_float), U3f(t_float), 
    U4f(t_float), U5f(t_float), U6f(t_float), 
    U7f(t_float), U8f(t_float), U9f(t_float), 
    U10f(t_float), U12f(t_float)
    
    #Assimilates function photosynthesis   
    A = Af(t)    
    
    # unpack y 
    T1, T2, C1, V1 = y[1], y[2], y[3], y[4]

     # compute RHS
    dy[1] = uconvert(u"K * s^-1" , rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12))
    dy[2] = uconvert(u"K * s^-1" , rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12))
    dy[3] = uconvert(u"mg * m^-3 * s^-1" ,rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10))
    dy[4] = uconvert(u"Pa* s^-1" ,rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7,U8, U9))

    #return [dy1, dy2, dy3, dy4]
end


    
