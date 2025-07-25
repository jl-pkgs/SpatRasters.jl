using SpatRasters, Dates, Test


@testset "sun_angle" begin
  lon = 120
  lat = 26.90248

  altitude, azimuth = angle_SunAzimuth(lat, DateTime(2010, 7, 25, 7, 47))
  @test altitude ≈ 31.976507857271315 # 30.65
  @test azimuth ≈ 82.80601450246094  # 82.08


  hours = 6:18
  N = length(hours)

  Hs = zeros(N)
  As = zeros(N)

  for i in 1:N
    time_local = DateTime(2025, 7, 25, hours[i])
    r = angle_SunAzimuth(lat, time_local)
    Hs[i], As[i] = r
  end
  As

  @test minimum(As) ≈ 72.44409880173262
  @test maximum(As) ≈ 287.5559011982674
  @test maximum(Hs) ≈ 82.63043358648318
end
