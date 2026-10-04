# In this version rhsPLL_PID! is updated to save all control values
# Probably this code will generate all data and save it as experimental 
# measurements. Th inference is made elsewhere. 

module HarverstSim

export harvestSimulation!
    # Required packages
    using CSV
    using DataFrames
    using Dates
    using Plots
    # using Unitful
    using DifferentialEquations
    using LaTeXStrings

    # weather data loader
    include("utils/weatherDataLoader.jl")
    # control functions Loader
    include("utils/controlLoader.jl")
    #psychrometrics
    include("models/psychrometrics.jl")
    # The RHS of greenhouse model
    include("utils/greenhouseModelGUI.jl")

    using .weatherDataLoader
    using .controls
    using .psychrometrics
    using .GreenhouseModel



    #---------------------------------------------------------------------------------------
    #
    # Mapping between simultion and global_dict
    row_map = Dict( "T1" => 1,"T2" => 2,"CO2" => 3,"VP" => 4, "A" => 7,"Tp" => 8,
        "Tout" => 11,"VPout" => 12, "PAR" => 13,"SW" => 14,"Wind" => 15,"Tsky" => 16,
        "Tsoil" => 17 
    )
    #---------------------------------------------------------------------------------------
    """

    One day simulation function

    """
    function  oneDaySimulation(p, u0, Weather, ctrls, evalTimes)
    
    t_span = (0, 24 * 3600.0)  
    # Define the ODE problem with current parameters `p`
    prob = ODEProblem((du, u, p, t) -> GreenhouseModel.simulationODErhs!(du, u, p, t, Weather, ctrls), u0, t_span, p)

    # Solve the ODE with time evaluation points
    sol = solve(prob, BS3(), saveat = evalTimes)

    dayunixGlobalTime = p.+ evalTimes 

    return sol, dayunixGlobalTime
    end

    #---------------------------------------------------------------------------------------
    # 
    """
    harvestSimulation function 
    """
    function harvestSimulation!(startDate, noDays, sim_dict)
        # auxiliary time variables  
        startUnixTime = Dates.datetime2unix(startDate)
        endDate = startDate + Dates.Day(noDays)
        daysRange = range(0, stop = noDays, length = noDays+1)
        sampleTime = 10 #minutes
        # sampling evaluation times  
        evalTimes = (0:60*sampleTime:24*3600)

        # weatherfile's path
        pathWeatherFile = "data/weatherData/weatherdata_NL_1988-07-11_2018-07.csv"
        # load weather interpolations
        Weather = weatherDataLoader.loadInterpolatedWeather(pathWeatherFile, startDate, endDate)

        # controlsorder's path 
        pathControlOrdersFile = "controlOrders/climateComputerOrders.json" 

        # initial conditions
        """
        format: 
        u0  = [ T1, T2, C1, V1, tempIntegralError, co2IntegralError, A, I3, U6, U8]
        units:[u"K", u"K", u"mg * m^-3", u"Pa", u"K s", u"mg * m^-3 * s", u"1", u"1", u"1", u"1" ]
        """
        u0 = [18+273.15, 23+273.15, 768.7, 1500.0, 0.0, 0.0, 0.0,0,0,0,0,0,0,0,0,0,0]

        
        for day_i in daysRange[1:  noDays]
            # Update starting daily simulation date 
            p = startUnixTime + 86400 * day_i  
            
            # load control interpolation
            controls.importControls(p, pathControlOrdersFile)
            ctrls =(
                ToutMax = controls.ToutMax, 
                Tset = controls.Tset,
                VentpBand = controls.VentpBand,
                ofset = controls.ofset,
                Light_on = controls.Light_on,
                CO2_set = controls.CO2_set)

            udia, tdia = oneDaySimulation(p, u0, Weather, ctrls, evalTimes)    

            u0 = udia[1:17,end]   

            
            allowed = Set(["T1", "T2", "Tp", "Tout", "Tsky", "Tsoil"])
            for (key, i) in row_map
                if key in allowed
                    append!(sim_dict[key], udia[i, :].-273.15)
                else
                    append!(sim_dict[key], udia[i, :])
                end 
            end
                append!(sim_dict["time"], tdia)
                append!(sim_dict["RH"], relative_humidity.(udia[row_map["T2"],:], udia[row_map["VP"],:]))
                append!(sim_dict["RHout"], relative_humidity.(udia[row_map["Tout"],:], udia[row_map["VPout"],:]))
        end    

        for key in keys(sim_dict)
            sim_dict[key][1] = sim_dict[key][2]
        end
        return sim_dict
    end


end
#---------------------------------------------------------------------------------------

# Main simulation routine
#---------------------------------------------------------------------------------------
"""
# Set up of simulations parameters
startDate = Dates.DateTime(2017, 9, 15, 0, 0)
noDays = 7

data_history = Dict(
    "time" => Float64[],
    "T1" => Float64[],
    "T2" => Float64[],
    "CO2" => Float64[],
    "VP" => Float64[],
    "A" => Float64[],
    "RH" => Float64[],
    "PAR" => Float64[],
    "Tp" => Float64[],
    "SW" => Float64[],
    "Tout" => Float64[],
    "Tsky" => Float64[],
    "Tsoil" => Float64[],
    "Wind" => Float64[],
    "RHout" => Float64[],
    "VPout" => Float64[]

)

data_history = harvestSimulation!(startDate, noDays, data_history)
"""