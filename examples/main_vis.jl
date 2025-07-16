using SpatRasters
using MakieLayers, GLMakie
import MakieLayers: imagesc, imagesc!

function imagesc(ra::SpatRaster, args...; kw...)
  lon, lat = st_dims(ra)
  imagesc(lon, lat, ra.A, args...; kw...)
end

function imagesc!(handle, ra::SpatRaster, args...; kw...)
  lon, lat = st_dims(ra)
  imagesc!(handle, lon, lat, ra.A, args...; kw...)
end

function plot_interp(fig, ra, X, Y; colors=amwg256, colorrange=(0, 60), title="IDW")
  ax, plt = imagesc!(fig, ra[:, :, 1]; colorrange, colors, title, force_show_legend=true)
  scatter!(ax, X[:, 1], X[:, 2]; color=Y, strokecolor=:white, strokewidth=0.6, colorrange, colormap=colors)
end
