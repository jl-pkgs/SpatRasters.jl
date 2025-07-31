using Test, SpatRasters
const sf = SpatRasters


@testset "earth_dist" begin
  p1 = sf.Point(100., 20.)
  p2 = sf.Point(101., 21.)

  r1 = earth_dist(p1, p2)
  r2 = earth_dist((100., 20.), (101., 21.))
  @test r1 == r2 == 152.5307854805719

  points = [
    100 20;
    101 20;
    100 21;
    101. 21.
  ]
  r = earth_dist((100., 20.), points)
  @test r[4] ≈ r2
end
