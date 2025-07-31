using SpatRasters, Test
const sf = SpatRasters

# p1 = sf.Point(1.0, 1.0)
# p2 = sf.Point(-1.5, 1.2)
# p3 = sf.Point(-1.2, -1.4)
# p4 = sf.Point(1.5, -1.3)

@testset "bilinear_irregular" begin
  p1 = sf.Point(-1., 1.)
  p2 = sf.Point(1., 2.0)
  p3 = sf.Point(-2., -1.0)
  p4 = sf.Point(2.0, -4.0)

  out_x, out_y = 0.0, 0.0
  t, s = get_fractional_distances_irregular(p1, p2, p3, p4, out_x, out_y)
  @test (t, s) == (0.375, 0.5)
end
