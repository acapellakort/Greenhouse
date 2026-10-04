# =============================================================================
# economics.jl  --  net-profit reward from AGC-2024 cost coefficients
# =============================================================================
# Turns a `simulate_twin` result (yield + metered resource use) into a euro
# net-profit breakdown. This is the Stage-3 reward signal. Coefficients from
# project note 08 (AGC-2024 simulator economics); cucumber price from the
# AGC-2018 Reference/Growers Production.csv (Prod_value_cum / Total_Prod_cum).

# --- price & cost coefficients (edit here) -----------------------------------
const PRICE_CUKE       = 0.889   # EUR / kg FW  (blended A+B, AGC-2018 Reference)
const COST_HEAT        = 0.09    # EUR / kWh    (natural gas heating)
const COST_CO2         = 0.30    # EUR / kg     (pure CO2, 300 EUR/ton)
const COST_ELEC_PEAK   = 0.30    # EUR / kWh    (07:00-23:00)
const COST_ELEC_OFF    = 0.20    # EUR / kWh    (23:00-07:00)
const FIXED_GH_PER_YR  = 15.00   # EUR / m2 / yr  greenhouse depreciation+maintenance
const PLANT_COST_STEM  = 0.22    # EUR / stem     (plant cost scales with density; ~2.5 stems -> 0.55/m2)

"""
    season_economics(res; price, cropdays, extra_fixed) -> NamedTuple

Net profit [EUR/m2] over the season, with a full breakdown. `res` is the
`simulate_twin` output. Revenue = Σ positive daily yield increments × price.
Variable costs = heating + CO2 + electricity (peak/off-peak lamps). Fixed costs
= greenhouse depreciation prorated by season length + plant cost (+ any extra).
"""
function season_economics(res; price = PRICE_CUKE,
                          cropdays = length(res.day), extra_fixed = 0.0)
    yld    = res.yield_FW
    dyield = length(yld) > 1 ? vcat(yld[1], diff(yld)) : copy(yld)
    revenue = sum(max.(dyield, 0.0)) * price

    heat_kWh = sum(res.heat_kWh)
    co2_kg   = sum(res.co2_kg)
    elec_pk  = sum(res.lamp_kWh_peak)
    elec_off = sum(res.lamp_kWh_off)

    heat_cost = heat_kWh * COST_HEAT
    co2_cost  = co2_kg   * COST_CO2
    elec_cost = elec_pk * COST_ELEC_PEAK + elec_off * COST_ELEC_OFF
    dens      = hasproperty(res, :density) && !isempty(res.density) ? maximum(res.density) : 2.5
    fixed     = (cropdays / 365.0) * FIXED_GH_PER_YR + PLANT_COST_STEM * dens + extra_fixed

    net = revenue - heat_cost - co2_cost - elec_cost - fixed
    return (; net, revenue, heat_cost, co2_cost, elec_cost, fixed,
             yield_kg = yld[end], heat_kWh, co2_kg, elec_kWh = elec_pk + elec_off)
end

"Pretty-print a season_economics breakdown."
function print_economics(e; label = "season")
    println("\n--- net profit ($label) ---")
    println("  revenue        ", round(e.revenue, digits = 2), "  EUR/m2   (", round(e.yield_kg, digits = 1), " kg/m2 @ ", PRICE_CUKE, " EUR/kg)")
    println("  - heating      ", round(e.heat_cost, digits = 2), "  EUR/m2   (", round(e.heat_kWh, digits = 1), " kWh/m2)")
    println("  - CO2          ", round(e.co2_cost, digits = 2), "  EUR/m2   (", round(e.co2_kg, digits = 1), " kg/m2)")
    println("  - electricity  ", round(e.elec_cost, digits = 2), "  EUR/m2   (", round(e.elec_kWh, digits = 1), " kWh/m2)")
    println("  - fixed        ", round(e.fixed, digits = 2), "  EUR/m2")
    println("  ------------------------------")
    println("  NET PROFIT     ", round(e.net, digits = 2), "  EUR/m2")
end
