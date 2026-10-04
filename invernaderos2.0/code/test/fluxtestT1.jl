using LaTeXStrings

# Variables climaticas
I1_t = LAIf(time_global)
I2_t = Iglobalf(time_global)
I4_t = Tskyf(time_global) 
I3_t = u[6,:] # Tuberia de calentamiento
I7_t = Tsoilf(time_global) 
I9_t = Idocelf(time_global)  

# Controles
U12_t = U12f(1.0).+ 0*time_global
U1_t = U1f(1.0).+ 0*time_global

# Variables de estado 
T1 = u[1,:]
T2 = u[2,:]
C1 = u[3,:]
V1 = u[4,:]


# Flujos 
r1_t = Rr1(I1_t,I2_t)
r5_t = Rr5(I2_t, U12_t)
r6_t =  Rr6(T1, I1_t, I3_t)
r7_t = Rr7(T1, I1_t, I4_t, U1_t)
h1_t = Rh1(T1, T2, I1_t)
l1_t = Rl1(T1, C1, V1, I1_t, I9_t)
r14_t = Rr14(T2, I7_t)
r15_t = Rr15(T1, I1_t, I4_t, U1_t) 


R1= r1_t .+ r5_t .+ r6_t .- h1_t .- l1_t .- r7_t.-  r14_t.- 0*r15_t

g1 = plot(r1_t, label=L"R_{PAR\_SunCan}", lw=2 )
plot!(r5_t, label=L"R_{NIR\_SunCan}", lw=2 )
plot!(r6_t, label=L"R_{PipeCan}", lw=2 )
plot!(-r7_t, label=L"R_{CanSky}", lw=2 )
plot!(-h1_t, label=L"H_{CanAir}", lw=2 )
plot!(-l1_t, label=L"L_{CanAir} ", lw=2 )
plot!(-r14_t, label=L"R_{CanFlr}", lw=2 )
#plot!(-r15_t, label=L"R_{CANCov,in}", lw=2 )
plot!(R1 , label=L"Total flux", lw =2)
display(g1)

g2= plot(T1.-273.15, label=L"T_{Can}", lw=2)
plot!(T2.-273.15, label=L"T_{air}", lw=2)
plot!(TTout.-273.15, label=L"T_{out}", lw=2, ylim=[0,50])
display(g2)

