# `temp_Brazil.csv`: https://github.com/rpkgs/spInterp
using SpatRasters, Statistics, Test
using RTableTools
include("main_vis.jl")

# indir = "$(@__DIR__)/../.." |> abspath
indir = "."
d = fread("$indir/data/prcp_st174_shiyan.csv")
b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
target = make_rast(; b, cellsize=0.01)

d = fread("$indir/data/temp_Brazil.csv")
b = bbox(-78, -34, -36, 5)
target = make_rast(; b, cellsize=0.25)

begin
  X = [d.lon d.lat] # d.alt
  Y = d[:, 3]
end

## 考虑角度与不考虑角度变化非常微弱
# Y = repeat(y, outer=(1, 24 * 30))
begin
  colorrange = (15, 30)
  radius = 500
  do_angle = false

  fig = Figure(; size=(1000, 1600))
  for (i, cdd) in enumerate([50, 100, 200, 400])
    for (j, nmax) in enumerate([4, 8])
      # @time ra_adw0 = interp(X, Y, target; wfun=weight_adw!, cdd, nmax, do_angle=false)
      @time ra = interp(X, Y, target; wfun=weight_adw!, cdd, nmax, radius, do_angle)
      plot_interp(fig[i, j], ra, X, Y; title="cdd = $cdd, nmax = $nmax", colorrange)
    end
  end
  save("Figure1_idw_(do_angle=$do_angle).png", fig)
  fig
end
