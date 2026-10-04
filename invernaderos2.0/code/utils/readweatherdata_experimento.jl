using CSV
using DataFrames
using Dates
using Unitful

# Function to compute Unix time and replace the column
function update_unix_time!(df)
    # Iterate over each row of the DataFrame
    for row in eachrow(df)
        # Construct a DateTime object from year, month, day, hour, and minute
        dt = DateTime(row.Year, row.Month, row.Day, row.Hour, row.Minute)
        
        # Convert the DateTime to Unix time
        unix_time = Dates.datetime2unix(DateTime(row.Year,row.Month,row.Day,row.Hour,row.Minute)) 
        
        # Substitute the Unix time into the first column
        row.Unix_Time = unix_time
    end
end

# Function to compute vapor preassure
function q2!(T)
    if T > 0
        return 611.21 * exp((18.678 - (T / 234.5)) * (T / (257.14 + T)))
    else
        return 611.21 * exp((23.036 - (T / 333.7)) * (T / (279.82 + T)))
    end
end

function fVaporPressure!(T,RH)
    return q2!(T)*(RH/100.0) 
end

function compute_Vapor_Pressure!(df)
    # Iterate over each row of the DataFrame
    for row in eachrow(df)
        # Construct a DateTime object from year, month, day, hour, and minute
        VP = fVaporPressure!(row.Temperature, row.Relative_Humidity)
        
        # Substitute the Unix time into the first column
        row.Vapor_Pressure = VP
    end
end


# Function to convert the CSV and filter between two given dates
function import_and_filter_csv(file_path, start_date, end_date)
    # Import the CSV file as a DataFrame
    df = CSV.read(file_path, DataFrame)

    # Rename the first column to "Unix_Time", and other columns for convenience
    rename!(df, 1  => "Unix_Time")
    rename!(df, 7  => "Temperature")
    rename!(df, 8  => "Relative_Humidity")
    rename!(df, 9  => "Vapor_Pressure")
    rename!(df, 13 => "External_Radiation")
    rename!(df, 14 => "Wind_Speed")
    
    # Update the value of Unix_Time conlumn 
    update_unix_time!(df)

    # Update the value of Vapor_Pressure conlumn 
    compute_Vapor_Pressure!(df)

    # Filter rows between the two given dates
    unix_start_date =  Dates.datetime2unix(start_date)
    unix_end_date =  Dates.datetime2unix(end_date)
    filtered_df =df[(df.Unix_Time .>= unix_start_date) .& (df.Unix_Time .<= unix_end_date), :]
    keep1 = ["Unix_Time", "Year", "Month", "Day", "Hour", "Minute", "Temperature", 
            "Relative_Humidity","Vapor_Pressure","External_Radiation","Wind_Speed"]
    filtered_df =select(filtered_df , keep1)

    filtered_df.Unix_Time = filtered_df.Unix_Time # u"s"
    filtered_df.Temperature = filtered_df.Temperature .+ 273.15 # u"K"
    filtered_df.Vapor_Pressure = filtered_df.Vapor_Pressure  # u"Pa"
    filtered_df.External_Radiation = filtered_df.External_Radiation #u"W * m^-2" 
    filtered_df.Wind_Speed = filtered_df.Wind_Speed .*  (5.0/18.0) #u"m * s^-1"

    return filtered_df
end

# Example usage:
file_path = "/Users/capella/Documents/trabajo/academia/unam/proyectos investigacion/bio/invernaderos/invernaderos2.0/code/data/dataset_meteo_experimento.csv"  # Path to your CSV file
 
filtered_data = import_and_filter_csv(file_path, start_date, end_date)

first(filtered_data,20)

# Create the functions to interpolate the weather data observations for all time
using DataInterpolations
using Plots

Toutf = DataInterpolations.ConstantInterpolation(filtered_data.Temperature, filtered_data.Unix_Time )
RHf   = DataInterpolations.ConstantInterpolation(filtered_data.Relative_Humidity, filtered_data.Unix_Time )
VPoutf   = DataInterpolations.ConstantInterpolation(filtered_data.Vapor_Pressure, filtered_data.Unix_Time )
Iglobalf  = DataInterpolations.ConstantInterpolation(filtered_data.External_Radiation, filtered_data.Unix_Time )
WindSpeedf = DataInterpolations.ConstantInterpolation(filtered_data.Wind_Speed, filtered_data.Unix_Time )

function CO2outf(t)
    return  834.7  #u"mg * m^-3"
end 

function LAIf(t)
    return 2 #u"m^2 * m^-2"
end

function Tpipef(f)
    return 15 + 273.15 #u"K"
end

function Tskyf(t)
    return -0.4+273.15 #"K"  
end
function TmechCoolf(t)
    return 273.15 #u"K"
end

function Tsoilf(t)
    return 18 + 273.15 #u"K"
end

Idocelf = Iglobalf

# I1 = LAI  
# I2 = Iglobal  
# I3 = Tpipef  +++++++?
# I4 = Tskif
# I5 = Tout
# I6 = TmechCool    
# I7 = Tsoil  
# I8 = WindSpeed m/s
# I9 = Idocel
#I10 = CO2out  
#I11 = VPout


# Temporalmente definimos Af la photosynthesis aqui
function Af(t)
    return 4 # u"mg * m^2 * s^-1"
end


# Use the interpolation function
# sample_time = Dates.datetime2unix(DateTime(1998, 7, 11, 9, 20))
# interpolated_value = Tout(sample_time)
# println("Interpolated value at $sample_time: $interpolated_value")
#plot(Tout)
#plot(RH)
#plot(VP)
#plot(Iglobal)
#plot(WindSpeed)
