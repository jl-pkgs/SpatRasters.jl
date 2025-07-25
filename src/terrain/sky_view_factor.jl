export SVF_azimuth, SVF_azimuth_simple, SVF


"slope in radian"
function dem_slope(p0::Point3{T}, p1::Point3{T}) where {T}
  dl = earth_dist(p0, p1) * 1000 # 水平面上的距离, [km] to [m]
  dz = p1.z - p0.z # [m]
  atan(dz / dl) # radians
end

function dem_slope(p0::Point3{T}, Points::Vector{Point3{T}}) where {T}
  map(p1 -> dem_slope(p0, p1), Points) # αs, H = pi/2 - maximum(αs)
end


function SVF_azimuth(Φ_sun::T, Φ_slope::T, β_slope::T, zenith_max::T) where {T<:Real}
  # Φ_sun = azimuth2deg(Φ_sun) # [天文学] -> [数学]，unnecessary, 二者统一即可
  # Φ_slope = azimuth2deg(Φ_slope) # [天文学] -> [数学]  
  Φ_sun = deg2rad(Φ_sun)
  Φ_slope = deg2rad(Φ_slope)

  H = zenith_max # radian
  sin(β_slope) * cos(Φ_sun - Φ_slope) * (H - sin(H) * cos(H)) + cos(β_slope) * sin(H)^2
end

function SVF_azimuth_simple(zenith_max::T) where {T<:Real}
  # H = pi / 2 - maximum(αs) # 天顶角 in [0, H]
  H = zenith_max # in radian
  sin(H)^2
end


"""
- `Φ_slope`: in [deg], 天文学方位角，正北为0
- `β_slope`: in [deg], 坡面角度
"""
function SVF(elev::SpatRaster, p0::Point; cellsize=nothing,
  radian=2.0, δψ=15, Φ_slope, β_slope, kernel::Function=SVF_azimuth)

  isnothing(cellsize) && (cellsize = st_cellsize(elev))

  β_slope = deg2rad(β_slope) # [deg] to [radian]


  z0 = st_extract(elev, [(p0.x, p0.y)]).value[1] # 
  P0 = Point3(p0.x, p0.y, z0)

  ∑ = 0.0
  n = 0 # 网格边界处，防止有的方向找不到网格
  ψs = δψ/2:δψ:360 # 天文学方位角
  N = length(ψs)
  for (i, Φ_sun) in enumerate(ψs)
    l = Line(; origin=p0, azimuth=Φ_sun, length=radian) # 200km^2
    points = intersect(elev, l; cellsize)
    length(points) == 0 && continue

    n += 1
    αs = dem_slope(P0, points) # 
    H = pi / 2 - maximum(αs) # 天顶角范围 [0, pi/2 - α]
    svf = kernel(Φ_sun, Φ_slope, β_slope, H)
    ∑ += svf
  end
  return ∑ / n # SVF, dΦ = 2pi/N
end


function SVF(elev::SpatRaster;
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
  rast(R, st_bbox(elev); bands=["SVF"])
end
