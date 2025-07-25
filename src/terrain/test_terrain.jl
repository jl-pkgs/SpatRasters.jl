using Dates, Test
using SpatRasters
using SpatRasters, ArchGDAL
using NetCDFTools

elev = rast("data/dem_etop01_G010deg.tif", FT=Float64)
@time MaxElevation = dem_angle_MaxElevation(ra; δψ=3, radian=2.0)

f = "./MaxElevation_etop01_G010deg.nc"
lon, lat = st_dims(f)
b = st_bbox(lon, lat)
@time MaxElevation = nc_read(f, "MaxElevation") # [radian]

ra_α = rast(MaxElevation, b)
# @profview R = dem_angle_MaxElevation(ra; δψ=3, radian=2.0)

time = DateTime(2010, 6, 10, 0)
@time R = sun_shade(ra_α, time)


using GLMakie, MakieLayers
using Shapefile

f_shp = "X:/rpkgs/sf.extract.R/inst/shp/Continents.shp"
shp = Shapefile.Table(f_shp)
imagesc(lon, lat, R.A)


begin
  time = DateTime(2010, 1, 10, 7)
  time_local = time + Hour(8)
  @show time_local

  @time R = sun_shade(ra_α, time)

  fig = Figure(; size=(1400, 800))
  
  limits = ((-180, 180), (-60, 90))
  # ticks = -8000:2000:8000 |> format_ticks
  ax, plt = imagesc!(fig[1, 1], lon, lat, R.A[:, :], axis=(; limits),
    # colorelevnge=(-8000, 8000), force_show_legend=true,
    # colorbar=(; ticks, width=20),
    title="(A) is_light")
  poly!(ax, shp.geometry; color=nan_color, strokecolor=:black, strokewidth=1.0)
  fig
end

begin
  elev = rad2deg.(MaxElevation[:, :, 60]) # 中午时刻最大仰角, [radian]
  fig = Figure(; size=(1400, 800))

  limits = ((-180, 180), (-60, 90))
  # ticks = -8000:2000:8000 |> format_ticks
  ax, plt = imagesc!(fig[1, 1], lon, lat, elev, axis=(; limits),
    # colorelevnge=(0, 8000), force_show_legend=true,
    # colorbar=(; ticks, width=20),
    title="(A) Altitude")
  poly!(ax, shp.geometry; color=nan_color, strokecolor=:black, strokewidth=1.0)
  fig
end

using SpatRasters: Point3
