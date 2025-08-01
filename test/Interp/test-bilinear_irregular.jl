using SpatRasters, Test
const sf = SpatRasters

# p1 = sf.Point(-1.5, 1.2)  # Upper left
# p2 = sf.Point(1.0, 1.0)   # Upepr right
# p3 = sf.Point(-1.2, -1.4) # Lower left
# p4 = sf.Point(1.5, -1.3)  # Lower right
pts_irregular = [
  Point(-1.5, 1.2),
  Point(1.0, 1.0),
  Point(-1.2, -1.4),
  Point(1.5, -1.3)
]

pts_vert_parallel = [
  Point(-1., 1.),
  Point(1., 2.),
  Point(-1., -1.),
  Point(1., -3.),
]

pts_hori_parallel = [
  Point(-1., 1.),
  Point(1., 3.),
  Point(-1., -1.),
  Point(1., -2.),
]

pts_both_parallel = [
  Point(-1., 1.),
  Point(1., 1.),
  Point(-1., -1.),
  Point(1., -1.),
]

@testset "bilinear FracDist" begin
  target = sf.Point(0, 0)
  out_x, out_y = 0.0, 0.0

  # 这样能解出来
  l_points = [
    "irregular" => pts_irregular,
    "vert" => pts_vert_parallel,
    "hori" => pts_hori_parallel,
    "both" => pts_both_parallel
  ]
  # r = frac_dist(pts_both_parallel..., out_x, out_y)
  @test frac_dist_parallellogram(pts_both_parallel[1:3]..., out_x, out_y) == (0.5, 0.5)
  @test frac_dist_parallellogram(pts_both_parallel[1:3]..., 0.2, 0.1) == (0.45, 0.6)

  for (k, points) in l_points[1:3]
    # printstyled(k, "\n", color=:green, bold=true)
    r1 = frac_dist_irregular(points..., out_x, out_y)
    r2 = frac_dist_uprights_parallel(points..., out_x, out_y)
    @test r1[1] ≈ r2[1]
    @test r1[2] ≈ r2[2]
  end

  for (k, points) in l_points
    r = frac_dist(points..., out_x, out_y)
    @test !isnan(r[1]) && !isnan(r[2])
  end
end
