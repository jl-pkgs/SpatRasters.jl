using Dates, Test
using SpatRasters
using SpatRasters, ArchGDAL
using NetCDFTools

elev = rast("data/dem_etop01_G010deg.tif", FT=Float64)
@time MaxElevation = dem_angle_MaxElevation(elev; δψ=3, radian=2.0)
# @profview R = dem_angle_MaxElevation(ra; δψ=3, radian=2.0)

lon, lat = st_dims(MaxElevation)
δψ = 3
ψs = δψ/2:δψ:360

A = Float32.(MaxElevation.A)
dims = (; lon, lat, ψs)
fout = "MaxElevation_etop01_G010deg_V2.nc"
ncsave(fout, true, (; units="radian");
  dims, MaxElevation=A)



using NaNStatistics
using NaNStatistics: nanmean

f = "MaxElevation_etop01_G010deg_V2.nc"
lon, lat = st_dims(f)
A = nc_read(f, "MaxElevation")
R = nanmean(A, dims=3)

# using GLMakie, MakieLayers
begin
  using CairoMakie, MakieLayers
  fig = Figure(; size=(1400, 800))
  imagesc!(fig, lon, lat, rad2deg.(R))
  save("Figure.png", fig)
end
