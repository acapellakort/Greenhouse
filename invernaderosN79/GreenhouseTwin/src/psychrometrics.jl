# =============================================================================
# psychrometrics.jl
# =============================================================================

"Saturation vapour pressure [Pa] at temperature T [K] (Buck-type formula over water / ice)."
function Pws(T)
    return ifelse(T > 273.15,
        611.21 * exp((18.678 - ((T - 273.15) / 234.5)) * ((T - 273.15) / (257.14 + T - 273.15))),
        611.15 * exp((23.036 - ((T - 273.15) / 333.7)) * ((T - 273.15) / (279.82 + T - 273.15))))
end

"Relative humidity [%] from air temperature [K] and vapour pressure [Pa]."
rhf(T, V) = 100 * V / Pws(T)
