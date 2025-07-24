using Dates
using SpatRasters
using GLMakie

begin
  lon = 120
  lat = 30.0

  hours = 0:24
  N = length(hours)

  Hs = zeros(n)
  As = zeros(n)

  for i in 1:N
    time_local = DateTime(2010, 6, 10, hours[i])
    r = angle_SunAzimuth(lat, time_local)
    Hs[i], As[i] = r
  end
end

points = map(azimuth, As)
scatter(points)
