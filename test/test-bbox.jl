@testset "bbox" begin
  lon = 70:140
  lat = 15:55
  b = st_bbox(lon, lat)
  @test b == bbox(69.5, 14.5, 140.5, 55.5)

  lon2, lat2 = bbox2dims(b; cellsize=1)
  @test length(lon2) == 71
  @test length(lat2) == 41
end
