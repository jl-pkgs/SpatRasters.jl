using SpatRasters, Statistics, Test
using RTableTools

# indir = "$(@__DIR__)/../.." |> abspath
indir = "."
d = fread("$indir/data/prcp_st174_shiyan.csv")

begin
  X = [d.lon d.lat] # d.alt
  Y = d.prcp
  # Y = repeat(y, outer=(1, 24 * 30))
  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)
end

neighbor = find_neighbor(target, X; radius=100, do_angle=true)
ra = interp(x, y, target)
# @profview r = interp_tps(x, Y, target; λ)

