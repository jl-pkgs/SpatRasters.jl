using Dates, Test, DataFrames
using SpatRasters

# using GLMakie
# DataFrame(; hour = hours, A = As)
# lon = 120
# lat = 26.90248

# t = DateTime(2015, 7, 25, 7, 47)
# time2local(t, 60, lon_ref=120.0) # 北京7点，lon=70° 3点
using SpatRasters, ArchGDAL


## 
ra = rast("data/dem_etop01_G010deg.tif", FT=Float64)
R = dem_angle_MaxElevation(ra; δψ=3, radian=2.0)
# altitude
