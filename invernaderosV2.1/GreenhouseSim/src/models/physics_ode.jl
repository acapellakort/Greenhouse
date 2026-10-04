"""
    rhsPLL_PID_sim!(du, u, p_struct, t)

The Right-Hand Side of the ODE system for the greenhouse simulation.
- `du`: Derivative vector
- `u`: State vector [T1, T2, C1, V1, integral_error]
- `p_struct`: Contains `.params` (ModelParams), `.controls_json` (ControlSet), 
              and `.controls_csv` (MeasuredControls)
- `t`: Time in seconds
"""
function rhsPLL_PID_sim!(du, u, p_struct, t)
    # 1. Unpack State
    T1, T2, C1, V1, integral_error = u
    
    # 2. Extract context
    params   = p_struct.params
    cs_json  = p_struct.controls_json
    cs_csv   = p_struct.controls_csv
    weather  = p_struct.weather # Assuming you store weather data here
    
    # 3. Unpack Parameters (via @unpack from Parameters.jl)
    # This keeps your math formulas clean and readable
    @unpack eta1, tau1, eta2 = params
    
    # 4. Handle Time
    t_calendar = t + params.start_time_unix
    current_date = Date(params.start_date + Second(round(Int, t)))

    # 5. Retrieve Weather/Environmental data
    # (Using the weather object we discussed earlier)
    I1  = weather[:LAI](t_calendar)
    I2  = weather[:Iglobal](t_calendar)
    I4  = weather[:Tsky](t_calendar)
    I5  = weather[:Tout](t_calendar)
    I6  = weather[:TmechCool](t_calendar)
    I7  = weather[:Tsoil](t_calendar)
    I8  = weather[:WindSpeed](t_calendar)
    I9  = (1 - eta1) * tau1 * eta2 * weather[:Idocel](t_calendar)
    I10 = weather[:CO2out](t_calendar)
    I11 = weather[:VPout](t_calendar)

    # 6. Photosynthesis assimilation calculation
    A = AcropFast(T2, C1, I9, I9) 

    # 7. Control variables
    # CSV-based controls (U1-U12, Tpipe)
    U1  = cs_csv.U1(t)
    U2  = cs_csv.U2(t)
    U3  = cs_csv.U3(t)
    U4  = cs_csv.U4(t)
    U5  = cs_csv.U5(t)
    U6  = cs_csv.U6(t)
    U7  = cs_csv.U7(t)
    U8  = cs_csv.U8(t)
    U9  = cs_csv.U9(t)
    U10 = cs_csv.U10(t)
    U11 = cs_csv.U11(t)
    U12 = cs_csv.U12(t)
    I3  = cs_csv.Tpipe(t)
    
    # JSON-based instructions (if needed for setpoints)
    # Tset = Controls.get_val(cs_json, :Tset, current_date, t)

    # 8. Compute ODEs
    # Using the rhs_fast function (must be in scope)
    du[1], du[2], du[3], du[4] = rhs_fast(T1, T2, C1, V1, I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, A, 
                                          U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U12)
    
    # 9. Update integral error and additional output
    du[5] = U11 - T2    # Integral of error for PI control
    u[5]  = integral_error   
    u[6]  = A           # Photosynthesis response
    
    return nothing
end