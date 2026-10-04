# Temporary definition of controls for code testing

function U1f(t) # pantalla termica
    return 0.0  
end

function U2f(t)  
    return 0.0  
end

function U3f(t) #Control del sistema de ventilador-almohadilla
    return 0.0  
end

function U4f(t) # Calentador de aire
    return 0.0  
end

function U5f(t) # pantalla sombreado
    return 0.0  
end

function U6f(t) 
    return 0.0  
end

function U7f(t) # vent forzada revisar leaking de calor nu4
    return 0.5  
end

function U8f(t) 
    return 0.5  
end

function U9f(t) 
    return 0.5 
end

function U10f(t) # Fuente de CO2
    return 0.0 #u"K"
end

function U11f(t) # Temp tuberia
    return 23 + 273.15 #u"K"
end

function U12f(t)
    return 1.0 #u"K"
end

