using SpatRasters, ArchGDAL

ra = rast("data/dem_etop01_G010deg.tif")
ra = rast(Float64.(ra.A), ra)

lon, lat = st_dims(ra)
# @profview 
R2 = SVF(ra; radian=2.0)

write_gdal(R2, "SVF_etop01_G010deg (radian=2.0deg).tif")

R5 = SVF(ra; radian=5.0)
write_gdal(R5, "SVF_etop01_G010deg (radian=5.0deg).tif")


# i, j = 1, 1
# p0 = Point(lon[i], lat[j])
# svf = SVF(ra, p0; radian=2.0, δψ=15, Φ_slope=0.0, β_slope=0.0, kernel=SVF_azimuth)

## 绘图 ------------------------------------------------------------------------
using MakieLayers, GLMakie
using Shapefile
GLMakie.activate!()

R2 = rast("SVF_etop01_G010deg (radian=2.0deg).tif")

f_shp = "D:/WSL/Ubuntu-20.04/rootfs/home/kong/R/x86_64-pc-linux-gnu-library/4.0/extract2/shp/Continents.shp"
shp = Shapefile.Table(f_shp)

# @profview
begin
  lon, lat = st_dims(R2)
  fig = Figure(; size=(850, 800))
  
  limits = ((-180, 180), (-60, 90))
  ticks = -8000:2000:8000 |> format_ticks
  
  ax, plt = imagesc!(fig[1, 1], lon, lat, ra.A[:, :], axis=(; limits),
    colorrange=(-8000, 8000), force_show_legend=true,
    colorbar=(; ticks, width=20),
    title="(A) DEM")
  poly!(ax, shp.geometry; color=nan_color, strokecolor=:black, strokewidth=1.0)

  ticks = 0.4:0.1:1.0 |> format_ticks
  ax, plt = imagesc!(fig[2, 1], lon, lat, R2.A[:, :], axis=(; limits),
    colorrange=(0.4, 1), force_show_legend=true, 
    colorbar=(; ticks, width=20), 
    title="(B) SVF")
  
  poly!(ax, shp.geometry; color=nan_color, strokecolor=:black, strokewidth=1.0)
  save("Figure1_SVF_(radian=2.0deg).jpg", fig)
  fig
end
