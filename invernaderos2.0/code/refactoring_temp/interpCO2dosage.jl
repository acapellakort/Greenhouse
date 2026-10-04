using Dates
using JSON
using DataInterpolations,  Plots
using CSV
using DataFrames


print(pwd())
file_path = "code/Data/CO2dosage.csv"
file_path_interp = "code/Data/CO2dosage_interp.csv"
df = CSV.read(file_path, DataFrame)
df_interp = CSV.read(file_path_interp, DataFrame)

print(pwd())
func = DataInterpolations.ConstantInterpolation(df.CO2_dosage, df.Times, extrapolate=true)
df_interp.CO2 = func(df_interp.time)

plot(df_interp.time,df_interp.CO2)

CSV.write(file_path_interp, df_interp)