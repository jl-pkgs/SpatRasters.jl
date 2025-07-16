using SpatRasters, Test

@testset "azimuth" begin
  @test angle_azimuth_sphere([0, 0], [0, 1]) == 0    # 正北，0
  @test angle_azimuth_sphere([0, 0], [1, 0]) == 90   # 正东，90
  @test angle_azimuth_sphere([0, 0], [0, -1]) == 180 # 正南，180
  @test angle_azimuth_sphere([0, 0], [-1, 0]) == 270 # 正西，270

  @test angle_azimuth([0, 0], [0, 1]) == 0    # 正北，0
  @test angle_azimuth([0, 0], [1, 0]) == 90   # 正东，90
  @test angle_azimuth([0, 0], [0, -1]) == 180 # 正南，180
  @test angle_azimuth([0, 0], [-1, 0]) == 270 # 正西，270
end


@testset "_weight_adw" begin
  dist = [10, 20, 40, 40, 50.]
  angle1 = [0, 10, 10, 180, 180.]
  angle0 = [0, 0, 0, 0, 0.]
  w0 = _weight_adw(dist, angle0)
  w1 = _weight_adw(dist, angle1)
  @test w1[4] > w1[3] # 对角获取更大的权重 
  @test sum(w0) ≈ 1
  @test sum(w1) ≈ 1
end
