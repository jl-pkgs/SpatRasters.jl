using SpatRasters, Statistics, Test
using RTableTools
include("main_vis.jl")

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

# neighbor = find_neighbor(target, X; radius=100, do_angle=true)
# weight_adw!(neighbor)
# neighbor.weight[:, :, 1]
# neighbor.angle[:, :, 1]
begin
  cdd = 100
  nmax = 5
  @time ra_adw0 = interp(X, Y, target; wfun=weight_adw!, cdd, nmax, do_angle=false)
  @time ra_adw1 = interp(X, Y, target; wfun=weight_adw!, cdd, nmax, do_angle=true)
  
  @time ra_tps1 = interp_tps(X, Y, target; λ = 0.01)
  @time ra_tps2 = interp_tps(X, Y, target; λ = 0.1)
  @time ra_idw = interp(X, Y, target; wfun=weight_idw!)
  # @profview ra_adw = interp(X, Y, target; wfun=weight_adw!)
  
  fig = Figure(; size=(1000, 800))
  plot_interp(fig[1, 1], ra_adw0, X, Y; title="ADW (angle=false)")
  plot_interp(fig[1, 2], ra_adw1, X, Y; title="ADW (angle=true)")
  plot_interp(fig[1, 3], ra_idw, X, Y; title="IDW")
  plot_interp(fig[2, 1], ra_tps1, X, Y; title="TPS (λ = 0.01)")
  plot_interp(fig[2, 2], ra_tps2, X, Y; title="TPS (λ = 0.1)")
  fig
end

@time ra_adw0 = interp(X, Y, target; wfun=weight_adw!, cdd=25, nmax=10, do_angle=false)
@time ra_adw1 = interp(X, Y, target; wfun=weight_adw!, cdd=25, nmax=10, do_angle=true)
