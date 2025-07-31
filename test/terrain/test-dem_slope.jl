using SpatRasters, Test, NaNStatistics, Dates
using Random

@testset "dem_angle_MaxElevation" begin
  Random.seed!(1234)

  b = bbox(100., 20., 102., 22.) # 5 * 5
  elev = make_rast(; b, cellsize=0.05)
  elev.A = rand(size(elev.A)...) * 10
  # display(elev.A[1:10, 1:10])

  # MaxElevation
  δψ = 3
  α = dem_angle_MaxElevation(elev; δψ)
  @test nanmaximum(α.A) <= 0.00188 # 0.0018745753805282214

  # SVF
  r = SVF(elev)
  @test minimum(r.A) >= 0.999

  # shade
  t = DateTime(2010, 6, 10, 11)
  l_shade = sun_shade(α, t)
  @test length(findall(.!l_shade.A)) == 40
end
