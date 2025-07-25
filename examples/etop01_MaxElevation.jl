using Dates, Test
using SpatRasters
using SpatRasters, ArchGDAL
using NetCDFTools

ra = rast("data/dem_etop01_G010deg.tif", FT=Float64)
@time MaxElevation = dem_angle_MaxElevation(ra; δψ=3, radian=2.0)
# @profview R = dem_angle_MaxElevation(ra; δψ=3, radian=2.0)

lon, lat = st_dims(MaxElevation)
δψ = 3
ψs = δψ/2:δψ:360

A = Float32.(MaxElevation.A)
dims = (; lon, lat, ψs)
ncsave("MaxElevation_etop01_G010deg.nc"; 
  dims, MaxElevation=A)



using NaNStatistics
using NaNStatistics: nanmean

f = "MaxElevation_etop01_G010deg.nc"
lon, lat = st_dims(f)
A = nc_read(f, "MaxElevation")
R = nanmean(A, dims=3)

using GLMakie, MakieLayers
imagesc(lon, lat, R)
