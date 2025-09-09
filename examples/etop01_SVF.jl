using SpatialRasterLite, ArchGDAL

elev = rast("data/dem_etop01_G010deg.tif", FT=Float64)
# @time MaxElevation = dem_angle_MaxElevation(elev; δψ=3, radian=2.0) # 每个方位角的最大仰角
# ψs = δψ/2:δψ:360
# write_gdal(elev, "MaxElevation_etop01_G010deg (vdian=2.0deg).tif")

# lon, lat = st_dims(elev)
R2 = SVF(elev; radian=2.0)
write_gdal(R2, "SVF_etop01_G010deg (radian=2.0deg)_V2.tif")

# R5 = SVF(elev; radian=5.0)
# write_gdal(R5, "SVF_etop01_G010deg (radian=5.0deg).tif")


# i, j = 1, 1
# p0 = Point(lon[i], lat[j])
# svf = SVF(elev, p0; radian=2.0, δψ=15, Φ_slope=0.0, β_slope=0.0, kernel=SVF_azimuth)

## 绘图 ------------------------------------------------------------------------
using MakieLayers, CairoMakie
using Shapefile
# using GLMakie
# GLMakie.activate!()

R2 = elevst("SVF_etop01_G010deg (radian=2.0deg)_V2.tif")

f_shp = "/mnt/x/rpkgs/sf.extract.R/inst/shp/Continents.shp"
shp = Shapefile.Table(f_shp)

function add_basemap!(ax)
  poly!(ax, shp.geometry; color=nan_color, strokecolor=:black, strokewidth=1.0)
end

# @profview
begin
  lon, lat = st_dims(R2)
  fig = Figure(; size=(850, 800))
  
  limits = ((-180, 180), (-60, 90))
  ticks = -8000:2000:8000 |> format_ticks
  
  ax, plt = imagesc!(fig[1, 1], lon, lat, elev.A[:, :], axis=(; limits),
    colorrange=(-8000, 8000), force_show_legend=true,
    colorbar=(; ticks, width=20),
    fun_axis=add_basemap!, title="(A) DEM")
  

  ticks = 0.4:0.1:1.0 |> format_ticks
  ax, plt = imagesc!(fig[2, 1], lon, lat, R2.A[:, :], axis=(; limits),
    colorrange=(0.99, 1), force_show_legend=true,
    colorbar=(; ticks, width=20), 
    fun_axis=add_basemap!, title="(B) SVF")  
  save("Figure1_SVF_(radian=2.0deg)_V2.png", fig)
  fig
end
