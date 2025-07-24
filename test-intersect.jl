# imagesc(lon, lat, R.A[:, :])
# i, j = 1, 1
# p0 = Point(lon[i], lat[j])
# svf = SVF(ra, p0; radian=2.0, δψ=15, Φ_slope=0.0, β_slope=0.0, kernel=SVF_azimuth)
# imagesc(lon, lat, A)

b = bbox(95.0, 15.0, 105.0, 25.0)
ra = make_rast(; b, cellsize=0.25)

p0 = Point(101.11, 20.11)

δψ = 15.0
ψs = δψ/2:δψ:360 # 天文学方位角
l_points = map(ψ -> begin
    l = Line(; origin=p0, azimuth=ψ, length=1.0) # 200km^2
    points = intersect(ra, l)
    @show ψ, length(points)
    @test length(points) >= 7
    # points
    (; points, ψ)
  end, ψs)



begin
  ψs = -360:15:360.
  p0 = Point(101.11, 20.11)
  
  l_points = map(ψ -> begin
      l = Line(; origin=p0, azimuth=ψ, length=2.0) # 200km^2
      points = intersect(ra, l)
      @show ψ, length(points)
      # @test length(points) >= 7
      # points
      (; points, ψ)
    end, ψs)
end


begin
  l = Line(; origin=Point(101.11, 20.11), azimuth=0.0)
  points = intersect(ra, l)
end

# l = Line(; origin=Point(101.11, 20.11), azimuth=0.0)
# points = intersect(ra, l)


## 判断交点
begin
  cols = resample_colors(amwg256, length(ψs))
  lon, lat = st_dims(ra)

  using GLMakie, MakieLayers
  fig = Figure(; size=(1000, 600))
  ax, plt = imagesc!(fig, lon, lat, ra.A)

  p0 = l.origin
  scatter!(ax, p0.x, p0.y; color=:black)

  for k = eachindex(l_points)
    (; points, ψ) = l_points[k]
    l = Line(; origin=Point(101.11, 20.11), azimuth=ψ, length=2.0) # 200km^2
    p0 = line_start(l)
    p1 = line_end(l)

    scatter!(ax, points; color=cols[k], markersize=3.5, strokecolor=cols[k])
    lines!([p0.x, p1.x], [p0.y, p1.y])
  end
  fig
  # save("Figure.png", fig)
end
