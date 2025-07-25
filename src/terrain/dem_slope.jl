export dem_angle_MaxElevation
export dem_slope


"slope in radian"
function dem_slope(p0::Point3{T}, p1::Point3{T}) where {T}
  dl = earth_dist(p0, p1) * 1000 # 水平面上的距离, [km] to [m]
  dz = p1.z - p0.z # [m]
  atan(dz / dl) # radians, [-90°, 90°]
end

function dem_slope(p0::Point3{T}, Points::Vector{Point3{T}}) where {T}
  map(p1 -> dem_slope(p0, p1), Points) # αs, H = pi/2 - maximum(αs)
end


## 提前算好，各个方向的最大坡度
function dem_angle_MaxElevation(elev::SpatRaster, p0::Point{T};
  δψ=3, radian=2.0, rastersize::RasterSize) where {T}

  z0 = st_extract(elev, [(p0.x, p0.y)]).value[1] # 
  P0 = Point3(p0.x, p0.y, z0)

  ψs = δψ/2:δψ:360 # 天文学方位角
  # ψs =[180.0]
  map(Φ_sun -> begin
      l = Line(; origin=p0, azimuth=Φ_sun, length=radian) # 200km^2
      points = intersect(elev, l, rastersize)
      length(points) == 0 && return NaN

      αs = dem_slope(P0, points) # 最大仰角对应的[radian]
      maximum(αs)
    end, ψs)
end


function dem_angle_MaxElevation(elev::SpatRaster; δψ=3, radian=2.0)
  # cellsize = st_cellsize(elev)
  rastersize = RasterSize(elev)
  lon, lat = st_dims(elev)
  nlon, nlat = length(lon), length(lat)
  # nlon = 20

  ψs = δψ/2:δψ:360 # 天文学方位角
  N = length(ψs)
  R = zeros(nlon, nlat, N)
  
  p = Progress(nlon)
  @inbounds @threads for i in 1:nlon
    next!(p)
    for j in 1:nlat
      p0 = Point(lon[i], lat[j])
      R[i, j, :] .= dem_angle_MaxElevation(elev, p0; δψ, radian, rastersize)
    end
  end
  rast(R, st_bbox(elev))
end
