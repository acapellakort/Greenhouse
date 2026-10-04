using CSV
using DataFrames
using Dates
using Unitful
using BenchmarkTools


# Include the script that imports constants
include("../utils/Version 2/import_constants.jl")

# include all functions
include("../utils/Version 2/climate_functionsV2.jl")

# Define reazonable values of all input variables 
T1 =  uconvert(u"K", 20u"°C")
T2 = uconvert(u"K", 22u"°C")
C1 = 680u"mg * m^-3"
V1 = 1200u"Pa"

A = 4u"mg * m^-2 * s^-1"

I1 = 2u"m^2 * m^-2"
I2 = 1000u"W * m^-2" 
I3 = uconvert(u"K", 15u"°C")
I4 = uconvert(u"K", -0.5u"°C")
I5 = uconvert(u"K", 23u"°C")
I6 = uconvert(u"K", 10u"°C")
I7 = uconvert(u"K", 17u"°C")
I8 = uconvert(u"m * s^-1", 23u"km * hr^-1")
I9 = 800u"W * m^-2" 
I10 = 810u"mg * m^-3"
I11 = 1500u"Pa"

U1  = 0.5
U2  = 0.5
U3  = 0.5
U4  = 0.5
U5  = 0.5
U6  = 0.5
U7  = 0.5
U8  = 0.5 
U9  = 0.5
U10 = 0.5
U11 = 0.5
U12 = 0.5

# ecuaciones para T1
r_1 = Rr1(I1, I2)
println("Value of r1: $r_1 ")
r_5 = Rr5(I2, U12)
println("Value of r5: $r_5 ")
r_6 = Rr6(T1, I1, I3)
println("Value of r6: $r_6 ")
r_7 = Rr7(T1, I1, I4, U1)
println("Value of r7: $r_7 ")
h_1 = Rh1(T1, T2, I1)
println("Value of h1: $h_1 ")
l_1 = Rl1(T1, C1, V1, I1, I9)
println("Value of l1: $l_1 ")

# ecuaciones para T2
h_1  = Rh1(T1, T2, I1) 
println("Value of h1: $h_1 ")
h_2  = Rh2(I5, U2) 
println("Value of h2: $h_2 ")
h_3  = Rh3(T2, V1, I6, U3)
println("Value of h3: $h_3 ")
h_4  = Rh4(T2, I3)
println("Value of h4: $h_4 ")
h_7  = Rh7(T2, I5, I8, U1, U5, U6, U7)
println("Value of h7: $h_7 ")
h_11 = Rh11(T2, I7)
println("Value of h11: $h_11 ")
h_12 = Rh12(U12)
println("Value of h12: $h_12 ")
r_8  = Rr8(I2)
println("Value of r8: $r_8 ")
l_2 = Rl2(U9)
println("Value of l2: $l_2 ")
r_10 = Rr10(T2, I1, I4, U1)
println("Value of r10: $r_10 ")


# Rhs de V1 
p_1  = Rp1(T1, C1, V1, I1, I9) 
println("Value of p1: $p_1 ")
p_2  = Rp2(U2) 
println("Value of p2: $p_2 ")
p_3  = Rp3(U9)
println("Value of p3: $p_3 ")
p_4  = Rp4(U4)
println("Value of p4: $p_4 ")
p_5  = Rp5(T2, V1, I5, I8, I11, U1, U5, U6, U7, U8)
println("Value of p5: $p_5 ")
p_6  = Rp6(T2, V1, U2)
println("Value of p6: $p_6 ")
p_7  = Rp7(T2, V1, I6, U3)
println("Value of p7: $p_7 ")


# Functiosn for C1
o_1  = Ro1(U4) 
println("Value of o1: $o_1 ") 
o_2  = Ro2(U10) 
println("Value of o2: $o_2 ") 
o_3  = Ro3(C1, I10, U2)
println("Value of o3: $o_3 ") 
o_4  = Ro4(A) 
println("Value of o4: $o_4 ") 
o_5  = Ro5(T2, C1, I5, I8, I10, U1, U5, U6, U7, U8)
println("Value of o5: $o_5 ") 



@btime rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12)
@btime rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)
@btime rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10)
@btime rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7, U8, U9)


rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12)
rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)
rhsC1(T2, C1, I5, I8, I10, A, U1, U2, U4, U5, U6, U7, U8, U10)
rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7, U8, U9)
