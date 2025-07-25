export sun_shade

function sun_shade(α::SpatRaster, time::DateTime; δψ=3)
  kmax = Int(360 / δψ)
  MaxElevation = α.A # 每个方位的最大仰角
  # cellsize = st_cellsize(elev)
  lon, lat = st_dims(α)
  nlon, nlat = length(lon), length(lat)

  R = zeros(nlon, nlat)
  # p = Progress(length(lon))
  @inbounds @threads for i in 1:nlon
    # next!(p)
    _lon = lon[i]
    for j in 1:nlat
      _lat = lat[j]
      
      time_local = time2local(time, _lon; lon_ref=0.0)
      H, A = angle_SunAzimuth(_lat, time_local; in_radian=true)
      H = max(H, 0) # 夜晚，设置H为0

      k = round(Int, rad2deg(A) / δψ) # 看落在哪一个方位角区间
      k = clamp(k, 1, kmax)

      _α = max(MaxElevation[i, j, k], 0)
      R[i, j] = H > _α
    end
  end
  rast(R, st_bbox(α); bands=["light"])
end

# ## 函数计算效率低，有内存优化空间
# function sun_shade(elev::SpatRaster, p0::Point{T}, time::DateTime;
#   radian=2.0, cellsize=nothing) where {T}
#   isnothing(cellsize) && (cellsize = st_cellsize(elev))

#   lon, lat = p0.x, p0.y
#   time_local = utc2local(time, lon)

#   z0 = st_extract(elev, [(p0.x, p0.y)]).value[1] # 
#   P0 = Point3(p0.x, p0.y, z0)

#   H, Φ_sun = angle_SunAzimuth(lat, time_local)

#   l = Line(; origin=p0, azimuth=Φ_sun, length=radian) # 200km^2
#   points = intersect(elev, l; cellsize)

#   αs = dem_slope(P0, points) # 判断最大s仰角
#   α = maximum(αs)
#   ## 太阳高度角>仰角，则光线能射入
#   H > α ? 1 : 0 # 1:亮; 0:阴
# end
