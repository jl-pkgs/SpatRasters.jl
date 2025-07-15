using SpatRasters, Statistics, Test
using RTableTools

indir = "$(@__DIR__)/../.." |> abspath
d = fread("$indir/data/prcp_st174_shiyan.csv")

@testset "tps" begin
  x = [d.lon d.lat d.alt]
  y = d.prcp
  Y = repeat(y, outer=(1, 2))

  λ = 0.01 # 小: 局部细节; 大: 空间上更加平缓
  @time tps1 = solve_tps(x, y, λ)
  @time tps2 = solve_tps(x, Y, λ)

  z1 = predict(tps1, x)
  z2 = predict(tps2, x)

  @test cor(z1, y)[1] .>= 0.999
  @test z2[:, 1] == z2[:, 2]
  @test maximum(abs.(z2[:, 1] - z1)) <= 1e-6
end
