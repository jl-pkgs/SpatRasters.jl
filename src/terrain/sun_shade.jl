export dem_angle_MaxElevation

## 提前算好，各个方向的最大坡度
function dem_angle_MaxElevation(elev::SpatRaster, p0::Point{T}; 
  δψ=3, radian=2.0, cellsize=nothing) where {T}
  
  isnothing(cellsize) && (cellsize = st_cellsize(elev))

  z0 = st_extract(elev, [(p0.x, p0.y)]).value[1] # 
  P0 = Point3(p0.x, p0.y, z0)

  ψs = δψ/2:δψ:360 # 天文学方位角
  map(Φ_sun -> begin
      l = Line(; origin=p0, azimuth=Φ_sun, length=radian) # 200km^2
      points = intersect(elev, l; cellsize)
      length(points) == 0 && return NaN

      αs = dem_slope(P0, points)
      maximum(αs)
    end, ψs)
end


function dem_angle_MaxElevation(elev::SpatRaster; δψ=3, radian=2.0)
  cellsize = st_cellsize(elev)
  lon, lat = st_dims(elev)
  nlon, nlat = length(lon), length(lat)

  ψs = δψ/2:δψ:360 # 天文学方位角
  N = length(ψs)
  R = zeros(nlon, nlat, N)

  p = Progress(length(lon))
  @threads for i in 1:nlon
    next!(p)
    for j in 1:nlat
      p0 = Point(lon[i], lat[j])
      R[i, j, :] = dem_angle_MaxElevation(elev, p0; δψ, radian, cellsize)
    end
  end
end


## TODO: 函数计算效率低，有内存优化空间
function sun_shade(elev::SpatRaster, p0::Point{T}, time::DateTime;
  radian=2.0, cellsize=nothing) where {T}
  isnothing(cellsize) && (cellsize = st_cellsize(elev))

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


function sun_shade(elev::SpatRaster, time::DateTime;)
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
        R[i, j] = sun_shade(elev, p0, time; cellsize, radia)
      catch ex
        @show i, j, ex
      end
    end
  end
  rast(R, st_bbox(elev); bands=["light"])
end


export sun_shade
