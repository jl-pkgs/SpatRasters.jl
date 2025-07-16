using SpatRasters, Test, RTableTools

@testset "interp" begin
  indir = "$(@__DIR__)/../.." |> abspath
  d = fread("$indir/data/prcp_st174_shiyan.csv")

  X = [d.lon d.lat] # d.alt
  Y = d.prcp
  # Y = repeat(y, outer=(1, 24 * 30))
  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)  
  # neighbor = find_neighbor(target, X; radius=100, do_angle=true)

  @time ra_idw = interp(X, Y, target; wfun=weight_idw!)
  @time ra_adw = interp(X, Y, target; wfun=weight_adw!, cdd=25, nmax=10, do_angle=true)
  @time ra_tps1 = interp_tps(X, Y, target; λ = 0.01)
  @time ra_tps2 = interp_tps(X, Y, target; λ = 0.1)
  # @profview ra_adw = interp(X, Y, target; wfun=weight_adw!)

  ## extract values
  points = st_points(X[:, 1], X[:, 2])
  z_idw = st_extract(ra_idw, points).value[:]
  z_adw = st_extract(ra_adw, points).value[:]
  z_tps1 = st_extract(ra_tps1, points).value[:]
  z_tps2 = st_extract(ra_tps2, points).value[:]

  @test cor(z_idw, z_adw) >= 0.9
  @test cor(z_idw, z_tps1) >= 0.9
  @test cor(z_idw, z_tps2) >= 0.9
end

# cor(Y, z_idw)
# cor(Y, z_adw)
# cor(Y, z_tps1)
# cor(Y, z_tps2)
