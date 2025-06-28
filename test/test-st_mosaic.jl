# using Test, Ipaper, Ipaper.sf, ArchGDAL

@testset "st_mosaic" begin
  bands = string.(1:4)
  r2 = rast(rand(4, 4, 4), bbox(-180.0, -30.0, 180.0, 0.0); bands)
  r1 = rast(rand(4, 4, 4), bbox(-180.0, -60.0, 180.0, -30.0); bands)

  rs = [r1, r2]
  r_big = st_mosaic(rs)
  @test st_bbox(r_big) == bbox(-180.0, -60.0, 180.0, 0.0)

  ## test for `merge_var`
  fs = ["r1.tif", "r2.tif"]
  write_gdal(r1, fs[1])
  write_gdal(r2, fs[2])

  box = st_bbox(st_bbox.(fs))
  R = merge_var(fs; box)
  @test cat(r2.A, r1.A, dims=2) ≈ R
  rm.(fs)
end
