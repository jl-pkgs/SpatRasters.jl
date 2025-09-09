using Dates, Test
using SpatialRasterLite
using SpatialRasterLite, ArchGDAL
using NetCDFTools

const sf = SpatialRasterLite
# elev = rast("data/dem_etop01_G010deg.tif", FT=Float64)
# @time MaxElevation = dem_angle_MaxElevation(elev; δψ=3, radian=2.0)
# p0 = sf.Point(90., 30.0)
# dem_angle_MaxElevation(elev, p0; δψ=3, radian=2.0)

function read_altitude(f)
  lon, lat = st_dims(f)
  b = st_bbox(lon, lat)
  MaxElevation = nc_read(f, "MaxElevation") # [radian]
  rast(MaxElevation, b)
end

f = "./MaxElevation_etop01_G010deg_V3.nc"
@time ra_α = read_altitude(f)
α_mean = rad2deg.(nanmean(ra_α.A, dims=3)[:, :, 1]) # 均值

time = DateTime(2010, 6, 10, 0)
@time R = sun_shade(ra_α, time)



using CairoMakie, MakieLayers
using Shapefile

# f_shp = "X:/rpkgs/sf.extract.R/inst/shp/Continents.shp"
f_shp = "/mnt/x/rpkgs/sf.extract.R/inst/shp/Continents.shp"
shp = Shapefile.Table(f_shp)

function add_basemap!(ax)
  poly!(ax, shp.geometry; color=nan_color, strokecolor=:black, strokewidth=1.0)
end

begin
  time = DateTime(2010, 1, 10, 0) + Minute(1.5*60)
  title = string(time)
  time_local = time + Hour(8)
  @show time_local
  @time R = sun_shade(ra_α, time)

  fig = Figure(; size=(1000, 1000))

  limits = ((-180, 180), (-60, 90))
  # ticks = -8000:2000:8000 |> format_ticks
  ax, plt = imagesc!(fig[1, 1], lon, lat, R.A[:, :], axis=(; limits),
    fun_axis=add_basemap!, title="(A) Day: $title")

  ax, plt = imagesc!(fig[2, 1], lon, lat, α_mean, axis=(; limits),
    colorrange=(0, 2), force_show_legend=true,
    fun_axis=add_basemap!, title="(B) Mean Altitude")

  # colorrange=(-8000, 8000), force_show_legend=true,
  # colorbar=(; ticks, width=20),
  save("Figure1_shade.png", fig)
  # fig
end
