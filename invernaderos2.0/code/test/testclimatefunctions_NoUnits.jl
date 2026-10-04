using CSV
using DataFrames
using Dates
using Unitful
using BenchmarkTools


# Include the script that imports constants
include("../utils/import_constants_NoUnits.jl")

# include all functions
include("../utils/climate_functionsV2_NoUnits.jl")
include("../utils/cropphoto_functions_NoUnits.jl")

# Define reazonable values of all input variables 
T1 = 20 + 273.15#u"K"
T2 = 22 + 273.15#u"K"
C1 = 680#u"mg * m^-3"
V1 = 1200#u"Pa"

#A = 4#u"mg * m^-2 * s^-1"

I1 = 2#u"m^2 * m^-2"
I2 = 1000#u"W * m^-2" 
I3 = 15 + 273.15#u"K"
I4 =-0.5 + 273.15#u"K"
I5 = 23 + 273.15#u"K"
I6 = 10 + 273.15#u"K"
I7 = 17 + 273.15#u"K"
I8 = 6.388#u"m * s^-1"
I9 = 800#u"W * m^-2" 
I10 = 610#u"mg * m^-3"
I11 = 1500#u"Pa"

Ac = Acrop(T1,C1,I9,I2)

U1  = 0.5
U2  = 0.5
U3  = 0.0 #Control del sistema de ventilador-almohadilla
U4  = 0.0
U5  = 0.5
U6  = 0.5
U7  = 0.5
U8  = 0.5 
U9  = 0.5
U10 = 0.0
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
o_4  = Ro4(Ac) 
println("Value of o4: $o_4 ") 
o_5  = Ro5(T2, C1, I5, I8, I10, U1, U5, U6, U7, U8)
println("Value of o5: $o_5 ") 

#s = rhs_fast(T1, T2, C1, V1, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, Ac, U1,  U2, U3,  U4, U5, U6, U7,U8, U9,U10, U12)
#print(s)

@btime rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12)
@btime rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)
@btime rhsC1(T2, C1, I5, I8, I10, Ac, U1, U2, U4, U5, U6, U7, U8, U10)
@btime rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7, U8, U9)

#println("Test the V1 function:")
#@btime      rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7, U8, U9)


T11=rhsT1(T1, T2, C1, V1, I1, I2, I3, I4, I9, U1, U12)
T21=rhsT2(T1, T2, I1, I2, I3, I4, I5, I6, I7, I8, U1, U2, U3, U5, U6, U7, U9, U12)
C11=rhsC1(T2, C1, I5, I8, I10, Ac, U1, U2, U4, U5, U6, U7, U8, U10)
V11=rhsV1(T1, T2, C1, V1, I1, I5, I6, I8, I9, I11, U1, U2, U3, U4, U5, U6, U7, U8, U9)

println("fast version")
@btime rhs_fast(T1, T2, C1, V1, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, Ac, U1,  U2, U3,  U4, U5, U6, U7,U8, U9,U10, U12)
println(s)

println([T11,T21,C11,V11])