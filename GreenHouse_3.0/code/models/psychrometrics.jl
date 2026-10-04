module psychrometrics

    export p_sat, relative_humidity, vapor_pressure, vapor_pressure_deficit

    """
        p_sat(T::Float64) ->  Float64

        Calculates the saturation water vapor pressure (in Pascals, Pa) based on the temperature `T` in Kelvin. 
        It uses different empirical formulas depending on whether the temperature is above or below 0°C (273.15K).

        # Inputs:
        - `T::Float64`: The temperature in Kelvin (K).

        # Outputs:
        - Returns the saturation water vapor pressure in Pascals (Pa).

        # Formula:
        - For T > 273.15K: 
            Pws = 611.21 * exp((18.678 - (T - 273.15) / 234.5) * ((T - 273.15) / (257.14 + T - 273.15)))
        - For T <= 273.15K: 
            Pws = 611.15 * exp((23.036 - (T - 273.15) / 333.7) * ((T - 273.15) / (279.82 + T - 273.15)))

        # References:
        - Stull, R. B. (2011). "Meteorology for Scientists and Engineers." Brooks/Cole, Cengage Learning.
    """       
    function p_sat(T::Float64)::Float64
        # Function to calculate the saturation water vapor pressure (in Pa) based on temperature T1.
        pws = ifelse.(T .> 273.15,
            611.21 .* exp.((18.678 .- ((T .- 273.15) ./ 234.5)) .* ((T .- 273.15) ./ (257.14 .+ T .- 273.15))),
            611.15 .* exp.((23.036 .- ((T .- 273.15) ./ 333.7)) .* ((T .- 273.15) ./ (279.82 .+ T .- 273.15)))
        )
        return pws
    end

    """
        relative_humidity(T::Float64, V::Float64) ->  Float64

    Calculates the relative humidity (RH) as a percentage based on the temperature `T` in Kelvin and the actual vapor pressure `V1` in Pascals. 
    The RH is computed by comparing the actual vapor pressure with the saturation vapor pressure at the given temperature.

    # Inputs:
    - `T::Float64`: The temperature in Kelvin (K).
    - `V::Float64`: The actual vapor pressure in Pascals (Pa).

    # Outputs:
    - Returns the relative humidity as a percentage (0 to 100%).

    # Formula:
    - RH = (V / p_sat(T)) * 100

    # References:
    - Lumsden, R. L. (2009). "Meteorology: An Introduction to the Atmosphere." Academic Press.
    """
    function relative_humidity(T::Float64, V::Float64)::Float64
        # Function to calculate relative humidity (%) using saturation water vapor pressure in Pa (Pws) and actual vapor pressure (V1) in Pa.
        return 100 .* V ./ p_sat(T)
    end

    """
        vapor_pressure_deficit(p_sat::Float64, rh::Float64) ->  Float64

    Calculates the vapor pressure deficit (VPD) in Pascals (Pa) based on the saturation vapor pressure (`VP_sat`) and the relative humidity (`rh`). 
    VPD represents the difference between the saturation vapor pressure and the actual vapor pressure at a given temperature.

    # Inputs:
    - `p_sat::Float64`: The saturation vapor pressure in Pascals (Pa).
    - `rh::Float64`: The relative humidity as a percentage (0 to 100%).

    # Outputs:
    - Returns the vapor pressure deficit (VPD) in Pascals (Pa).

    # Formula:
    - vapor_pressure_deficit = p_sat * (1 - rh / 100)

    # References:
    - Jones, H. G. (1992). "Plants and Microclimate: A Quantitative Approach to Environmental Plant Physiology." Cambridge University Press.
    """
    function vapor_pressure_deficit(p_sat::Float64, rh::Float64)::Float64
        # Function to calculate vapor pressure deficit (VPD) in Pa.
        return p_sat .* (1 .- rh ./ 100)
    end

    """
        vapor_pressure(T::Float64, rh::Float64) ->  Float64

    Calculates the vapor pressure (VP) in Pascals (Pa) based on the temperature (T)
    # Inputs:
    - `T::Float64`: Temperature (K).
    - `rh::Float64`: The relative humidity as a percentage (0 to 100%).

    # Outputs:
    - Returns the vapor pressure in Pascals (Pa).

    # Formula:
    - vapor_pressure = p_sat(T) * (rh / 100)

    # References:
    - Jones, H. G. (1992). "Plants and Microclimate: A Quantitative Approach to Environmental Plant Physiology." Cambridge University Press.
    """
    function vapor_pressure(T::Float64,rh::Float64)::Float64
        return p_sat(T)*(rh/100.0) 
    end


end 


"""
using .psychrometrics

psat = p_sat(300.0)

println(psat)
"""