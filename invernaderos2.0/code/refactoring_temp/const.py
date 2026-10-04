#%%
def Struct( typ, varid, prn, desc, units, val, ok='1'):
    return { 'val':val,'units':str(units), 'latex':prn, 'desc':desc, 'comment':ok }


ALPHA ={
    ################## alpha ##################
    'alpha1': Struct(typ='Cnts', varid='alpha1', prn= '$\alpha_1$',
                    desc="Heat capacity of one square meter of the canopy", units="J * K**-1 * m**-2", val=1.5e3, ok='Se regreso al valor original'), # Capacidad calórifica de un m^2 de dosel (theta[0])
    'alpha2': Struct(typ='Cnts', varid='alpha2', prn= '$\alpha_2$',
                    desc="Global NIR absorption coefficient of the canopy", units="1", val=0.35, ok="ok"), # Coeficiente global de absorción NIR del dosel
    'alpha3': Struct(typ='Cnts', varid='alpha3', prn= '$\alpha_3$',
                    desc="Surface of the heating pipe", units="m**2*m**-2", val=0.3), # Superficiedelatuber ́ıadecalentamiento
    'alpha4': Struct(typ='Cnts', varid='alpha4', prn= '$\alpha_4$',
                    desc="Convection heat exchange coefficient of canopy leaf to greenhouse air", units="W * m**-2 * K**-1", val=15, ok="ok"), # Coeficiente de intercambio de calor por conveccio ́n de la hoja del dosel al aire del invernadero
    'alpha5': Struct(typ='Cnts', varid='alpha5', prn= '$\alpha_5$',
                    desc="Specific heat capacity of greenhouse air", units="J * K**-1 * kg**-1", val=1.5e3, ok="ok"), # Capacidad calorıfica especifica delaire del invernadero
    'alpha6': Struct(typ='Cnts', varid='alpha6', prn= '$\alpha_6$',
                    desc="Greenhouse floor surface area", units="m**2", val=1e4, ok="ok"), # Área de la superficie del piso del invernadero
    'alpha7': Struct(typ='Cnts', varid='alpha7', prn= '$\alpha_7$',
                    desc="Global NIR absorption coefficient of the floor", units="1", val=0.5, ok="ok"), # Coeficiente global de absorcio ́n NIR del piso
    'alpha8': Struct(typ='Cnts', varid='alpha8', prn= '$\alpha_8$',
                    desc="PAR absorption coefficient of the cover", units="1", val=1, ok='no dan el valor'), # Coeficiente de absorci ́on PAR de la cubierta # En el artículo no dan el valor
    'alpha9': Struct(typ='Cnts', varid='alpha9', prn= '$\alpha_9$',
                    desc="NIR absorption coefficient of the cover", units="1", val=1, ok='no dan el valor'),
    'alpha12': Struct(typ='Cnts', varid='alpha12', prn= '$\alpha_{12}$',
                    desc="Total lamp radiation per square meter ", units = "W * m**-2", val=74, ok='Calculo Antonio')                
}


BETA = {
    ################## beta ##################
    'beta1': Struct(typ='Cnts', varid='beta1', prn= '$\beta_1$',
                    desc="Canopy extinction coefficient for PAR radiation", units="1", val=0.7, ok="ok"),
    'beta2': Struct(typ='Cnts', varid='beta2', prn= '$\beta_2$',
                    desc="Extinction coefficient for PAR radiation reflected from the floor to the canopy", units="1", val=0.7, ok="ok"), 
    'beta3': Struct(typ='Cnts', varid='beta3', prn= '$\beta_3$',
                    desc="Canopy extinction coefficient for NIR radiation", units="1", val=0.27, ok="ok")
}


GAMMA = {
    ################## gamma ##################
    'gamma':  Struct(typ='Cnts', varid='gamma', prn= '$\gamma$',
                    desc="Psychometric constan", units="Pa * K**-1", val=65.8, ok="ok"),  # Constante psicrom ́etrica #ok 
    'gamma1': Struct(typ='Cnts', varid='gamma1', prn= '$\gamma_1$',
                    desc="Length of the heating pipe", units="m * m**-2", val=1.25, ok='ok, usé el valor de Texas'), # Longitud de la tubería de calentamiento (Almería)
    'gamma2': Struct(typ='Cnts', varid='gamma2', prn= '$\gamma_2$',
                    desc="Latent heat of water evaporation", units="J * kg_water**-1", val=2.45e6, ok="ok"), # Calor latente de evaporaci ́on del agua #ok
    'gamma3': Struct(typ='Cnts', varid='gamma3', prn= '$\gamma_3$',
                    desc="Strength of boundary layer of canopy for vapor transport", units="s * m**-1", val=275, ok="ok"), # Resistencia de la capa l ́ımite del dosel para transporte de vapor # ok
    'gamma4': Struct(typ='Cnts', varid='gamma4', prn= '$\gamma_4$',
                    desc="Minimum stomatal resistance of the canopy", units="s * m**-1", val=82.0, ok="ok"), # Resistenciaestom ́aticam ́ınimadeldosel # ok
    'gamma5': Struct(typ='Cnts', varid='gamma5', prn= '$\gamma_5$',
                    desc="Slope of the differentiable switch for the stomatal resistance model", units="m * W**-2", val=-1, ok="ok")
}


DELTA = {
    ################## delta ##################
    'delta1': Struct(typ='Cnts', varid='delta1', prn= '$\delta_1$',
                    desc="Radiation above the canopy that defines sunrise and sunset", units="W * m**-2", val=5, ok="ok"), # Radiación por encima del dosel que define el amanecer y la puesta de sol # ok
    'delta2': Struct(typ='Cnts', varid='delta2', prn= '$\delta_2$',
                    desc="Empirically determined parameter", units="W * m**-2", val=4.3, ok="ok"), # Parámetro determinado empíricamente # ok
    'delta3': Struct(typ='Cnts', varid='delta3', prn= '$\delta_3$',
                    desc="Empirically determined parameter", units="W * m**-2", val=0.54, ok="ok"), # Parámetro determinado empíricamente # ok
    'delta4': Struct(typ='Cnts', varid='delta4', prn= '$\delta_4$',
                    desc="Coefficient of the CO2 transpiration in the day", units="ppm**-2", val=6.1e-7, ok="ok"), 
    'delta5': Struct(typ='Cnts', varid='delta5', prn= '$\delta_5$',
                    desc="Coefficient of the CO2 transpiration in the night", units="ppm**-2", val=1.1e-11, ok="ok"), 
    'delta6': Struct(typ='Cnts', varid='delta6', prn= '$\delta_6$',
                    desc="Coefficient of the vapour pressure in the day", units="Pa**-2", val=4.3e-6, ok="ok"),   
    'delta7': Struct(typ='Cnts', varid='delta7', prn= '$\delta_7$',
                    desc="Coefficient of the vapour pressure in the night", units="Pa**-2", val=5.2e-6, ok="ok")
}


EPSIL = {
    ################## epsilon ##################
    'epsil1': Struct(typ='Cnts', varid='epsil1', prn= '$\epsilon_1$',
                    desc="FIR emission coefficient of the heating pipe", units="1", val=0.88, ok="ok"), # Coeficiente de emisión FIR de la tubería de calentamiento # ok 
    'epsil2': Struct(typ='Cnts', varid='epsil2', prn= '$\epsilon_2$',
                    desc="Canopy FIR emission coefficient", units="1", val=1, ok="ok"), # Coeficiente de emisión FIR del dosel # ok
    'epsil3': Struct(typ='Cnts', varid='epsil3', prn= '$\epsilon_3$',
                    desc="Sky FIR emission coefficient", units="1", val=1, ok="ok"), # Coeficiente de emisión FIR del cielo # ok
    'epsil4': Struct(typ='Cnts', varid='epsil4', prn= '$\epsilon_4$',
                    desc="Floor FIR emission coefficient", units="1", val=1, ok="ok"), # Coeficiente de emisión FIR del piso # ok
    'epsil5': Struct(typ='Cnts', varid='epsil5', prn= '$\epsilon_5$',
                    desc="Thermal screen FIR emission coefficient", units="1", val=1, ok='?'), # Coeficiente de emisión FIR de la pantalla térmica
    'epsil6': Struct(typ='Cnts', varid='epsil6', prn= '$\epsilon_6$',
                    desc="External cover FIR emission coefficient", units="1", val=0.44, ok='ok,usé el valor de Texas')
}


ETA = {
    ################## eta ##################
    'eta1':  Struct(typ='Cnts', varid='eta1', prn= '$\eta_1$',
                    desc="Proportion of global radiation that is absorbed by greenhouse building elements", units="1", val=0.1, ok="ok"),  # Proporción de la radiación global que es absorbida por los elementos de construcción del invernadero # ok
    'eta2':  Struct(typ='Cnts', varid='eta2', prn= '$\eta_2$',
                    desc="Ratio between PAR radiation and external global radiation", units="1", val=0.5, ok="ok"),  # Razón entre la radiación PAR y la radiación global externa ¿0.5?
    'eta3':  Struct(typ='Cnts', varid='eta3', prn= '$\eta_3$',
                    desc="Ratio between NIR radiation and global external radiation", units="1", val=0.5, ok="ok"),  # Razón entre la radiación NIR y la radiación global externa # ok 
    'eta4':  Struct(typ='Cnts', varid='eta4', prn= '$\eta_4$',
                    desc="Conversion factor for CO2 of mg*m**−3 to ppm", units="ppm * mg**-1 * m**3", val=0.554, ok="ok"),  # Factor de conversión de mg m−3 CO2 a ppm # ok 
    'eta5':  Struct(typ='Cnts', varid='eta5', prn= '$\eta_5$',
                    desc="Fan-pad system efficiency", units="1", val=0.5, ok="ok"),  # Eficiencia del sistema de ventilador-almohadilla # no da el valor en el articulo
    'eta6':  Struct(typ='Cnts', varid='eta6', prn= '$\eta_6$',
                    desc="Ventilation power reduction factor", units="m**3 * m**-2 * s**-1", val=1, ok='Falta valor'),  # Factor de reduccio ́n de la potencia de ventilación # Falta valor
    'eta7':  Struct(typ='Cnts', varid='eta7', prn= '$\eta_7$',
                    desc="Ratio between ceiling ventilation area and total ventilation area", units="1", val=0.5, ok='no dan valor en el artículo'),  # Razón entre el área de ventilación en el techo y el  área de ventilación total  # no da el valor en el articulo
    'eta8':  Struct(typ='Cnts', varid='eta8', prn= '$\eta_8$',
                    desc="Ratio between ceiling and total ventilation area, if there is no chimney effect", units="1", val=0.9, ok="ok"),  # Razón entre el área de ventilación techo y total, si no hay efecto de chimenea # ok
    'eta9':  Struct(typ='Cnts', varid='eta8', prn= '$\eta_9$',
                    desc="", units="1", val=0, ok='No esta en el código'),  # Razón entre el área de ventilación lateral y el área de ventilación total # no hay eta9
    'eta10': Struct(typ='Cnts', varid='eta10', prn= '$\eta_{10}$',
                    desc="Shadow effect on the discharge coefficient", units="1", val=0, ok='Falta valor, en los ejemplos del artículo no se considera'), # Efecto de la sombra sobre el coeficiente de descarga # Falta valor
    'eta11': Struct(typ='Cnts', varid='eta11', prn= '$\eta_{11}$',
                    desc="Effect of shadow on the global wind pressure coefficient", units="1", val=0, ok='falta valor'), # Efecto de la sombra sobre el coeficiente de presión global del viento # Falta valor, aunque en los ejemplos del artículo no se considera
    'eta12': Struct(typ='Cnts', varid='eta12', prn= '$\eta_{12}$',
                    desc="Amount of vapor that is released when a joule of sensible energy is produced by the direct air heater", units="kg_vapour * J**-1", val=4.43e-8, ok="ok"), # Cantidad de vapor que es liberado cuando un joule de energía sensible es producido por el calentador de aire directo # ok
    'eta13': Struct(typ='Cnts', varid='eta13', prn= '$\eta_{13}$',
                    desc="Amount of CO2 that is released when a joule of sensible energy is produced by the direct air heater", units="mg_CO2 * J**-1", val=0.057, ok="ok"),
    'eta14': Struct(typ='Cnts', varid='eta14', prn= '$\eta_{14}$',
                    desc="Percentage lamps radiation that is NIR", units="1", val=0.18, ok='Aaron dio el valor '),
    'eta15': Struct(typ='Cnts', varid='eta15', prn= '$\eta_{15}$',
                    desc="Porcentaje de la radiacion de las lamparas que es calor en el reflector", units="1", val=0.0074, ok='Aaron dio el valor '),
    'eta16': Struct(typ='Cnts', varid='eta16', prn= '$\eta_{16}$',
                    desc="Porcentaje de la radiacion de las lamparas que es calor directo al aire del invernadero", units="1", val=0.39, ok='Aaron dio el valor '),
    'eta17': Struct(typ='Cnts', varid='eta17', prn= '$\eta_{17}$',
                    desc="Percentage lamps radiation that is PAR", units="1", val=0.3626, ok='Aaron dio el valor ')
    
}

LAMB = {
    ################## lamb ##################
    'lamb1': Struct(typ='Cnts', varid='lamb1', prn= '$\lambda_1$',
                    desc="Performance coefficient of the mechanical acceleration system", units="1", val=0, ok='Falta valor, en los ejemplos del artículo no se considera'), # Coeficiente de desempen ̃o del sistema de enfriamiento meca ́nico # Falta valor, aunque en los ejemplos del artículo no se considera
    'lamb2': Struct(typ='Cnts', varid='lamb2', prn= '$\lambda_2$',
                    desc="Electrical capacity of the mechanical cooling system", units="W", val=0, ok='Falta valor, en los ejemplos del artículo no se considera'), # Capacidad el ́ectrica del sistema de enfriamiento meca ́nico # Falta valor, aunque en los ejemplos del artículo no se considera
    'lamb3': Struct(typ='Cnts', varid='lamb3', prn= '$\lambda_3$',
                    desc="Convictive heat exchange coefficient between soil and greenhouse air", units="W * m**-2 * K**-1", val=1, ok='1 . Falta valor, en los ejemplos del artículo no se considera'), # Coeficiente de intercambio de calor convictivo entre el suelo y el aire del invernadero # Falta valor, aunque en los ejemplos del artículo no se considera
    'lamb4': Struct(typ='Cnts', varid='lamb4', prn= '$\lambda_4$',
                    desc="Heat capacity of direct air heater", units="W", val=5*(10**5), ok='Dr Antonio dio el valor'), # Capacidad calor ́ıfica del calentador de aire directo
    'lamb5': Struct(typ='Cnts', varid='lamb5', prn= '$\lambda_5$',
                    desc="Cover surface", units="m**2", val=1.8e4, ok='ok,tomé el valor de Holanda, el de Texas es muy grande (9e4)'), # Superficie de la cubierta # ok --> tomé el valor de Holanda, el de Texas es muy grande (9e4)
    'lamb6': Struct(typ='Cnts', varid='lamb6', prn= '$\lambda_6$',
                    desc="Variable of heat exchange by convection between the roof and the outside air", units="W * m_cover**-2 * K**-1", val=2.8, ok=' 2.8 ok, usé el valor de Texas'), # Variable de intercambio de calor por convecci ́on entre la cubierta y el aire exterior # ok ---> usé el valor de Texas
    'lamb7': Struct(typ='Cnts', varid='lamb7', prn= '$\lambda_7$',
                    desc="Variable of heat exchange by convection between the roof and the outside air", units="J * m**-3 * K**-1", val=1, ok=' 1 ok, usé el valor de Texas'), # Variable de intercambio de calor por convecci ́on entre la cubierta y el aire exterior # ok ---> usé el valor de Texas
    'lamb8': Struct(typ='Cnts', varid='lamb8', prn= '$\lambda_8$',
                    desc="Variable of heat exchange by convection between the roof and the outside air", units="1", val=2, ok='ok,usé el valor de Texas'),
    'lamb9':Struct(typ='Cnts',varid='lamb9',prn= '$\lambda_9$', 
                    desc='eficiencia del intercambio de calor entre las capas de aire', units="1", val=0.05, ok="ok")
}


RHO = {
    ################## rho ##################
    'rho1':Struct(typ='Cnts', varid='rho1', prn= '$\rho_1$',
                    desc="PAR reflection coefficient", units="1", val=0.07,ok="ok"),
    'rho2': Struct(typ='Cnts', varid='rho2', prn= '$\rho_2$',
                    desc="Floor reflection coefficient PAR", units="1", val=0.65,ok = "ok"), 
    'rho3': Struct(typ='Cnts', varid='rho3', prn= '$\rho_3$',
                    desc="Air density", units="kg * m**-3", val= 1.2,ok = 'El valor es el de la densidad del aire al nivel del mar'),
    'rho4': Struct(typ='Cnts', varid='rho4', prn= '$\rho_4$',
                    desc="Air density see level", units="kg * m**-3", val= 1.2,ok = 'El valor es el de la densidad del aire al nivel del mar')
}


TAU = {
    'tau1': Struct(typ='Cnts', varid='tau1', prn= '$\tau_1$',
                    desc="PAR transmission coefficient of the Cover", units="1", val=1,ok = 'En el artículo no dan su valor'),
    'tau2': Struct(typ='Cnts', varid='tau2', prn= '$\tau_2$',
                    desc="FIR transmission coefficient of the Cover", units="1", val=1, ok ='En el artículo no dan su valor'),
    'tau3': Struct(typ='Cnts', varid='tau3', prn= '$\tau_3$',
                    desc="FIR transmission coefficient of the thermal screen", units="1", val=0.11,ok = 'ok --> usé el valor de Texas')
}


NU ={
    'nu1': Struct(typ='Cnts', varid='nu1', prn= '$\nu_1$',
                    desc="Shadowless discharge coefficient", units="1", val=0.65,ok = "ok"), 
    'nu2': Struct(typ='Cnts', varid='nu2', prn= '$\nu_2$',
                    desc="Global wind pressure coefficient without shadow", units="1", val=0.1,ok="ok"),
    'nu3': Struct(typ='Cnts', varid='nu3', prn= '$\nu_3$',
                    desc="Side surface of the greenhouse", units="m**2", val=3,ok  = 'En ejemplos del artículo usan valor cero,estaba en 900'), 
    'nu4': Struct(typ='Cnts', varid='nu4', prn= '$\nu_4$',
                    desc="Leakage coefficien", units="1", val=1e-4,ok="ok"), 
    'nu5': Struct(typ='Cnts', varid='nu5', prn= '$\nu_5$',
                    desc="Maximum ceiling ventilation area", units="m**2", val=2e3, ok = ' 0.2*alpha6 --> ok'), 
    'nu6': Struct(typ='Cnts', varid='nu6', prn= '$\nu_6$',
                    desc="Vertical dimension of a single open respirator", units="m", val=1,ok="ok"), 
    'nu7': Struct(typ='Cnts', varid='nu7', prn= '$\nu_7$',
                    desc="Soil thermal conductivity", units="W * m**-1 * K**-1", val=0.85, ok="ok"), 
    'nu8': Struct(typ='Cnts', varid='nu8', prn= '$\nu_8$',
                    desc="Floor to ground distance", units="m", val=0.64,ok="ok")
}


PHI = {
    ################## phi ##################
    'phi1': Struct(typ='Cnts', varid='phi1', prn= '$\phi_1$',
                    desc="External diameter of the heating pipe", units="m", val=51e-3,ok = "ok"),
    'phi2': Struct(typ='Cnts', varid='phi2', prn= '$\phi_2$',
                    desc="Average height of greenhouse air", units="m", val=4, ok = 'Se regreso a valor original'), 
    'phi3': Struct(typ='Cnts', varid='phi3', prn= '$\phi_3$',
                desc='Masa molar del aire', units="kg * kmol**-1", val=28.9647, ok = "ok" ), 
    'phi4': Struct(typ='Cnts', varid='phi4', prn= '$\phi_4$',
                desc='Altitud del invernadero', units="metros", val=0, ok = "ok"),
    'phi5': Struct(typ='Cnts', varid='phi5', prn= '$\phi_5$',
                    desc="Water vapor contained in the fan-pad system", units="kg_water * kg_air**-1", val=0.014, ok="ok"), 
    'phi6': Struct(typ='Cnts', varid='phi6', prn= '$\phi_6$',
                    desc="Water vapor contained in the outside air", units="kg_water * kg_air**-1", val=0.0079, ok='este es el valor correcto a 21 grados C y 50 % de HR'), # Vapor de agua contenido en el aire exterior
    'phi7': Struct(typ='Cnts', varid='phi7', prn= '$\phi_7$',
                    desc="Capacity of air flow through the pad", units="m**3 * s**-1", val=16.7,ok = "ok"), 
    'phi8': Struct(typ='Cnts', varid='phi8', prn= '$\phi_8$',
                    desc="Air flow capacity of forced ventilation system", units="m**3 * s**-1", val=666.6, ok = 'https://farm-energy.extension.org/greenhouse-ventilation/'),
    'phi9': Struct(typ='Cnts', varid='phi9', prn= '$\phi_9$',
                    desc="Fog system capacity", units="kg * s**-1", val=0.916, ok='Como en Holanda')
}


PSI = {
    ################## psi ##################
    'psi1':Struct(typ='Cnts', varid='psi1', prn= '$\psi_1$',
                    desc="Molar mass of water", units="kg * kmol**-1", val=18,ok = "ok"), 
    'psi2': Struct(typ='Cnts', varid='psi2', prn= '$\psi_2$',
                    desc="Capacity of the external CO2 source", units="mg * s**-1", val=1.33*(10**4),ok="ok"),
    'psi3': Struct(typ='Cnts', varid='psi3', prn= '$\psi_3$',
                    desc="Molar mass of the CH2O", units="g * mol_CH2O**-1", val=30.031,ok = "ok")
}


OMEGA = {
    ################## omega ##################
    'omega1': Struct(typ='Cnts', varid='omega1', prn= '$\omega_1$',
                    desc="Gravity acceleration constant", units="m * s**-2", val=9.81, ok = "ok"), 
    'omega2': Struct(typ='Cnts', varid='omega2', prn= '$\omega_2$',
                    desc="Molar gas constant", units="J * kmol**-1 * K**-1", val= 8.314e3, ok = "ok"), 
    'omega3': Struct(typ='Cnts', varid='omega3', prn= '$\omega_3$',\
                        desc="Percentage of CO2 absorbed by the canopy", units=" 1 ", val=0.03/4.0,ok = 'Deberia depender del modelo de la planta')
}

# %%
