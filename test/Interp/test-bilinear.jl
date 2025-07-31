using SpatRasters, Test


@testset "bilinear" begin
  lon = 70:5:140
  lat = 55:-5:15
  b = st_bbox(lon, lat)
  z = rand(length(lon), length(lat), 2)
  ra_low = rast(z, b)

  Lon = 70:2.5:140
  Lat = 55:-2.5:15
  B = st_bbox(Lon, Lat)
  Z = bilinear(lon, lat, z, Lon, Lat; na_rm=true)
  ra_high = rast(Z, B)
  @test size(Z) == (29, 17, 2)

  ra_bi = interp(ra_low, ra_high; method="bilinear")
  ra_nr = interp(ra_low, ra_high; method="nearest")
  @test ra_bi.A == Z
  @test size(ra_bi) == size(ra_nr)
end
