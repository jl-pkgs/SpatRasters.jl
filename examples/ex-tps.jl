using SpatRasters, Statistics, Test
using RTableTools
using Distances
using NearestNeighbors

# indir = "$(@__DIR__)/../.." |> abspath
begin
  indir = "."
  d = fread("$indir/data/prcp_st174_shiyan.csv")
  x = [d.lon d.lat] # d.alt
  y = d.prcp
  X = collect(x')

  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)
  domain = st_coords(target)
end


begin
  Y = repeat(y, outer=(1, 24*30))
  λ = 0.01 # 小: 局部细节; 大: 空间上更加平缓
  @time r = interp_tps(x, Y, target; λ)
end
# @profview r = interp_tps(x, Y, target; λ)

begin
  limits = bbox2lims(b)
  lon, lat = st_dims(target)

  using GLMakie, MakieLayers
  fig = Figure(; size=(800, 600))

  ax, plt = imagesc!(fig[1, 1], lon, lat, r.A[:, :, 1];
    colors=amwg256, colorrange=(0, 60),
    title="λ=$λ", # limits,
    force_show_legend=true, col_rev=false, axis=(; limits))

  sc = scatter!(ax, x[:, 1], x[:, 2], color=y;
    strokecolor=:white, strokewidth=0.8, 
    colormap=amwg256, colorrange=(0, 60))
  fig
end
