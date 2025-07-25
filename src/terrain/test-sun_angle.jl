using Dates, Test, DataFrames
using SpatRasters
using GLMakie

DataFrame(; hour = hours, A = As)
lon = 120
lat = 26.90248

t = DateTime(2015, 7, 25, 7, 47)
time2local(t, 60, lon_ref=120.0) # 北京7点，lon=70° 3点
