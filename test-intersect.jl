using SpatRasters, Test
using ArchGDAL
using MakieLayers
using Ipaper

ra = rast("data/dem_etop01_G010deg.tif")
ra = rast(Float64.(ra.A), ra)

lon, lat = st_dims(ra)
nlon, nlat = length(lon), length(lat)

begin
  A = zeros(size(ra.A))
  p = Progress(length(lat))
  
  for i in 1:nlon
    next!(p)
    for j in 1:nlat
      p0 = Point(lon[i], lat[j])

      try
        svf = SVF(ra, p0;
          radian=2.0, δψ=15, Φ_slope=0.0, β_slope=0.0, kernel=SVF_azimuth)
        A[i, j] = svf

      catch ex
        @show i, j, ex
      end
    end
  end
end

i, j = 1, 1

imagesc(lon, lat, A)



p0 = Point(101.11, 20.11)

  
l_points = map(ψ -> begin
    l = Line(; origin=p0, azimuth=ψ, length=1.0) # 200km^2
    points = intersect(ra, l)
    @show ψ, length(points)
    @test length(points) >= 7
    # points
    (; points, ψ)
  end, ψs)

begin
  # b = bbox(95.0, 15.0, 105.0, 25.0)
  # ra = make_rast(; b, cellsize=0.25)

  ψs = -360:15:360.
  p0 = Point(101.11, 20.11)
  l_points = map(ψ -> begin
      l = Line(; origin=p0, azimuth=ψ, length=2.0) # 200km^2
      points = intersect(ra, l)
      @show ψ, length(points)
      @test length(points) >= 7
      # points
      (; points, ψ)
    end, ψs)
end

# begin
#   l = Line(; origin=Point(101.11, 20.11), azimuth=ψ)
#   points = intersect(ra, l)
# end
# l = Line(; origin=Point(101.11, 20.11), azimuth=0.0)
# points = intersect(ra, l)
cols = resample_colors(amwg256, length(ψs))

begin
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
