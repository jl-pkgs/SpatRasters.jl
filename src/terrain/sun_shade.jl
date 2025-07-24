## TODO: 函数计算效率低，有内存优化空间
function sun_shade(elev::SpatRaster, p0::Point{T}, time::DateTime; radian=2.0) where {T}
  cellsize = st_cellsize(elev)

  lon, lat = p0.x, p0.y
  time_local = utc2local(time, lon)

  z0 = st_extract(elev, [(p0.x, p0.y)]).value[1] # 
  P0 = Point3(p0.x, p0.y, z0)

  H, Φ_sun = angle_SunAzimuth(lat, time_local)

  l = Line(; origin=p0, azimuth=Φ_sun, length=radian) # 200km^2
  points = intersect(elev, l; cellsize)

  αs = dem_slope(P0, points) # 判断最大s仰角
  α = maximum(αs)
  ## 太阳高度角>仰角，则光线能射入
  H > α ? 1 : 0 # 1:亮; 0:阴
end


function sun_shade(elev::SpatRaster;
  radian=2.0, δψ=15, Φ_slope=0.0, β_slope=0.0, kernel=SVF_azimuth)

  cellsize = st_cellsize(elev)
  lon, lat = st_dims(elev)
  nlon, nlat = length(lon), length(lat)

  R = zeros(size(elev.A))
  p = Progress(length(lon))
  @threads for i in 1:nlon
    next!(p)
    for j in 1:nlat
      p0 = Point(lon[i], lat[j])
      try
        R[i, j] = SVF(elev, p0; cellsize, radian, δψ, Φ_slope, β_slope, kernel)
      catch ex
        @show i, j, ex
      end
    end
  end
  rast(R, st_bbox(elev); bands=["light"])
end


export sun_shade
