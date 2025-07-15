using SpatRasters, Statistics, Test
using RTableTools

# indir = "$(@__DIR__)/../.." |> abspath
indir = "."
d = fread("$indir/data/prcp_st174_shiyan.csv")

begin
  x = [d.lon d.lat] # d.alt
  X = collect(x')
  y = d.prcp
  Y = repeat(y, outer=(1, 24 * 30))

  ## 空间插值
  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)
  # λ = 0.01 # 小: 局部细节; 大: 空间上更加平缓
  # @time r = interp_tps(x, Y, target; λ)
end

neighbor = find_neighbor(target, x; radius=100)
ra = interp(x, y, target)
# @profview r = interp_tps(x, Y, target; λ)

