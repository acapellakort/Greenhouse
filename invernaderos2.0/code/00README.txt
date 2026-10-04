Codigo del invernadero:

Implementado:

-- Modelo de clima simplificado
     a) El modelo de clima opera por día, al 
     final del día se actuializan lños controles.
     Con este diseño el modelo de crecimiento puede 
     actuar una vez al dia, en el loop exterior. 
     b) El modelo de fotosintesis actua instantanemament
     c) Método de solución del sistema de EDO:
     Se usa Tsit5(), este método no es el mas rápido, 
     pero el sistema parace tener problemas para converger 
     si se usan otros. 

    # Trapezoid()  8.571 ms (333725 allocations: 5.19 MiB)   <- No converge
    # QNDF()      23.523 ms (819563 allocations: 12.65 MiB)  <- No converge
    # BS3()       37.827 ms (1499374 allocations: 23.06 MiB) 33 s 28.886 s 
    # FBDF()      48.172 ms (1689795 allocations: 26.02 MiB) <- No converge
    # Tsit5()     53.546 ms (2105611 allocations: 32.38 MiB) 55 s 
    # DP5()       61.941 ms (2325378 allocations: 35.75 MiB)
    # Vern7()     82.261 ms (3332109 allocations: 51.30 MiB)

-- Modelo de fotosintesis
    Se programo con resistencia estomática constante. 
    FALTA: modelo con dependencia en variables ambientales.

Desempeño del modelo:

Tiempo:

-- Modelo completo 90 dias:               55 s
-- Modelo con fotosintesis constante:     33 s
    
Obs: El desempeño del modelo se hizo más lento con las funciones
de control por horas.  

-- El modelo de fotosintesis toma 2.745μs este tiempo 
es muy lento ya que como se ejecuta en serie incrementa
el tiempo de ejecución del RHS para cada evaluación.  

Con la mejora en las funciones, ahora el tiempo el desempeño es:
1.082 μs (61 allocations: 976 bytes)
y el del modelo completo es:
28.886 s (887603334 allocations: 13.33 GiB)
19.6 s con la implementacion del fast rhs este se evalua en 5.2 ms 

-- La evaluacion del RHS del modelo de clima es:


 ToDo: - Implementar en un única función todo el modelo de 
       fotosintesis. Esto ayudara con los tiempos totales. 
       - Mejorar el desarrollo de los rhs de T1, T2, C1, V1
  T1    1.739 μs (91 allocations: 1.42 KiB)
  T2    1.854 μs (98 allocations: 1.53 KiB)
  C1    0.995 μs (61 allocations: 976 bytes)
  V1    2.439 μs (128 allocations: 2.00 KiB)


  5 segs 10 segs



  Parametros a inferir:


    "alpha1": {
        "val": 1500.0,
        "units": "J * K^-1 * m^-2",
        "desc": "Heat capacity of one square meter of the canopy",
        "latex": "$\u0007lpha_1$",
        "comment": "Se regreso al valor original"
    },
    "alpha2": {
        "val": 0.35,
        "units": "1",
        "desc": "Global NIR absorption coefficient of the canopy",
        "latex": "$\u0007lpha_2$",
        "comment": "ok"
    },       

  "alpha4": {                   <-------------------- este es el primero
        "val": 5,
        "units": "W * m^-2 * K^-1",
        "desc": "Convection heat exchange coefficient of canopy leaf to greenhouse air",
        "latex": "$\u0007lpha_4$",
        "comment": "estaba 15, perp en Vanthor es 5"}


    "beta2": {
        "val": 0.7,
        "units": "1",
        "desc": "Extinction coefficient for PAR radiation reflected from the floor to the canopy",
        "latex": "$\beta_2$",
        "comment": "ok"
    },

        "gamma3": {
        "val": 275,
        "units": "s * m^-1",
        "desc": "Strength of boundary layer of canopy for vapor transport",
        "latex": "$\\gamma_3$",
        "comment": "ok"
    },

    "gamma4": {
        "val": 82.0,
        "units": "s * m^-1",
        "desc": "Minimum stomatal resistance of the canopy",
        "latex": "$\\gamma_4$",
        "comment": "ok"
    },

 "nu4": {
        "val": 1e-2,
        "units": "1",
        "desc": "Leakage coefficien",
        "latex": "$\nu_4$",
        "comment": "ok valor 1e-4 en Vanthor lo hacemos mas grande para estabilizar el sistema de ODE"
    },

        "nu7": {
        "val": 0.85,
        "units": "W * m^-1 * K^-1",
        "desc": "Soil thermal conductivity",
        "latex": "$\nu_7$",
        "comment": "ok"
    },

    g_eff() effective stomatal resistence (funcion)

(1 - eta1) * tau1   como uan so0la constante 


    "eta1": {
        "val": 0.1,
        "units": "1",
        "desc": "Proportion of global radiation that is absorbed by greenhouse building elements",
        "latex": "$\\eta_1$",
        "comment": "ok"
    },



    "tau1": {
        "val": 0.77,
        "units": "1",
        "desc": "PAR transmission coefficient of the Cover",
        "latex": "$\tau_1$",
        "comment": "En el art\u00edculo no dan su valor"
    },



Definir una fotosinthesis efectiva:

def A2(C,I, alpha_R=0.0625, beta_R=-4, alpha_f=0.0622,beta_a=30):
    '''
    Units 
    [ A ] = mu_mol m**-2 s**-1
    [ I ] = mu_mol m**-2 s**-1
    [ C ] = mu_bar = ppm

    Values:
    alpha_R =  (20+5)/400 = 0.0625 (Fig 7 (b) ajuste lineal)
    beta_R = -4 mu_mol m**-2 s**-1 (Fig 1) 
    alpha_f = 0.0622 (Fig 8)
    beta_a = 30 mu_mol m**-2 s**-1 (Fig 7 (a) and (b) )

    Source: Yin, X. & Struik, P.C.,  C3 and C4 photosynthesis 
    models: An overview from the perspective of crop modelling, 
    Wageningen Journal of Life Sciences, 2009/12/01, 
    https://doi.org/10.1016/j.njas.2009.07.001 

    '''
    return min([alpha_R*C+beta_R,alpha_f*I,beta_a])

    inferir los valores., Este experimento es facil. 
    Se4 simula con el modelo completo y se hace la inferencia de regreso. 
    Luego se compara 



plot(rh_data)


function ExpU(alpha)
    tcan_p, tair_p, rh_p, co2_p = forwardMap!(alpha)
    logL = -sum((rh_data.-rh_p).^2)./sigma[3].^2-sum((tair_data.-tair_p).^2)/sigma[3]^2
    return logL
end

function ExpSupp(alpha)
    if 1 .< alpha &&  alpha.< 25 
        return true
    else
        return false
    end
end

Exp = jtwalk( n=1, U=ExpU, Supp=ExpSupp)

Run!(Exp, T=500, x0=2.0, xp0=11.0)

ToDo:
Set up a data adquisition routine to save all control variables.
These will be used as control inputs into the model to set up the inverse problem. 

3OCT24

Parametro: alpha4

La inferencia de este paramtro no es clara. Al correr el MCMC con disttintos
niveles de ruido la cadena oscila



gamma4: funciona!! y muy bien

