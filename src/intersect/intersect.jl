export Point, Line, intersect
export line_start, line_end
export azimuth2slope
export SVF_azimuth, SVF_azimuth_simple
export SVF
import Base: intersect

abstract type AbstractPoint{T} end

Base.@kwdef mutable struct Point3{T} <: AbstractPoint{T}
  x::T
  y::T
  z::T
end

Base.@kwdef mutable struct Point{T} <: AbstractPoint{T}
  x::T
  y::T
end

"""
方位角转换为数学角度，限定在0-360°

```bash
0 -> 90
90 -> 0
180 -> -90
270 -> -180
```
"""
function azimuth2deg(x)
  x = -x + 90
  (x > 360) && (x -= 360)
  (x < 0) && (x += 360)
  x
end

azimuth2slope(ψ) = tan(deg2rad(azimuth2deg(ψ)))

Base.@kwdef mutable struct Line{T}
  origin::Point{T} = Point(0.0, 0.0)
  azimuth::T = 0.0 # deg, 正北为0, 顺时针为正
  length::T = 2.0 # in deg
  k::T = azimuth2slope(azimuth)
end


@inline line_start(line::Line) = line.origin

function line_end(line::Line)
  p0 = line.origin
  length = line.length
  x, y = p0.x, p0.y
  θ = line.azimuth |> azimuth2deg |> deg2rad
  Point(x + cos(θ) * length, y + sin(θ) * length)
end

function st_bbox(line::Line)
  p0 = line_start(line)
  p1 = line_end(line)
  xmax = max(p0.x, p1.x)
  xmin = min(p0.x, p1.x)
  ymax = max(p0.y, p1.y)
  ymin = min(p0.y, p1.y)
  bbox(; xmin, ymin, xmax, ymax)
end

function is_vertical(line::Line; eps=1e-4)
  abs(mod(line.azimuth, 180)) <= eps # 0或180，认为是垂线
end

"和x轴的交点，如果k = Inf or -Inf，则需要用y轴相交的方法"
function intersect_x(line::Line{T}, xs::Vector{T}) where {T}
  (; k) = line
  p0 = line_start(line)
  x0, y0 = p0.x, p0.y
  map(x -> begin
      y = k * (x - x0) + y0
      Point(x, y)
    end, xs)
end

function intersect_y(line::Line{T}, ys::Vector{T}) where {T}
  (; k) = line
  p0 = line_start(line)
  x0, y0 = p0.x, p0.y
  map(y -> begin
      x = (y - y0) / k + x0
      Point(x, y)
    end, ys)
end

function intersect(ra::SpatRaster, line::Line)
  # lon, lat = st_dims(ra)
  cellx, celly = st_cellsize(ra)
  b = st_bbox(ra)
  lon = b.xmin:cellx:b.xmax # 采用的是网格边界
  lat = celly < 0 ? (b.ymax:celly:b.ymin) : (b.ymin:celly:b.ymax)

  bl = st_bbox(line)
  ## 判断交点
  ilat = bl.ymin .<= lat .<= bl.ymax # 
  ilon = bl.xmin .<= lon .<= bl.xmax # 判断大致的范围
  xs = lon[ilon]
  ys = lat[ilat]

  points_y = intersect_y(line, ys) # 与水平线的交点
  points_x = intersect_x(line, xs) # 与垂线的交点

  points = cat(points_y, points_x, dims=1) |> rm_empty

  x = map(p -> p.x, points)
  inds = sortperm(x)
  points = @view points[inds] # 对points进行排序  

  ## 然后两点判断一个网格位置
  _cellij(ra, points)
end


function earth_dist(p1::Point3{T}, p2::Point3{T}) where {T}
  earth_dist((p1.x, p1.y), (p2.x, p2.y))
end

"slope in radian"
function cal_α(p0::Point3{T}, p1::Point3{T}) where {T}
  dl = earth_dist(p0, p1) * 1000 # 水平面上的距离, [km] to [m]
  dz = p1.z - p0.z # [m]
  atan(dz / dl) # radians
end


function cal_α(p0::Point3{T}, Points::Vector{Point3{T}}) where {T}
  map(p1 -> cal_α(p0, p1), Points) # αs, H = pi/2 - maximum(αs)
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
function SVF(ra::SpatRaster, p0::Point;
  radian=2.0, δψ=15, Φ_slope, β_slope, kernel::Function=SVF_azimuth)

  β_slope = deg2rad(β_slope) # [deg] to [radian]

  ψs = δψ/2:δψ:360 # 天文学方位角
  N = length(ψs)

  z0 = st_extract(ra, [(p0.x, p0.y)]).value[1] # 
  P0 = Point3(p0.x, p0.y, z0)

  ∑ = 0.0
  n = 0 # 网格边界处，防止有的方向找不到网格
  for (i, Φ_sun) in enumerate(ψs)
    l = Line(; origin=p0, azimuth=Φ_sun, length=radian) # 200km^2
    points = intersect(ra, l)
    length(points) == 0 && continue
    # @show i, Φ_sun, length(Points)
    # @show points, length(points)
    n += 1
    αs = cal_α(P0, points) # 
    H = pi / 2 - maximum(αs) # 天顶角范围 [0, pi/2 - α]
    svf = kernel(Φ_sun, Φ_slope, β_slope, H)
    ∑ += svf
  end
  return ∑ / n # SVF, dΦ = 2pi/N
end


"""
仅用于计算相交

# Arguments
points: 与网格边界相交的所有点

# Return
每两个点 确定一个网格中心，返回的是网格中心的[x, y, elev]
"""
function _cellij(ra::SpatRaster, points::AbstractVector{Point{T}}) where {T}
  n = length(points)
  b = st_bbox(ra)
  lon, lat = st_dims(ra)
  cellx, celly = st_cellsize(ra)
  nx, ny = size(ra)[1:2]

  map(i -> begin
      p1 = points[i]
      p2 = points[i+1]
      x = (p1.x + p2.x) / 2
      y = (p1.y + p2.y) / 2
      p = _location_fast((x, y); b, cellx, celly, nx, ny) # (i, j)

      if isnothing(p)
        nothing # 没有找到点
      else
        i, j = p
        Point3(lon[i], lat[j], ra.A[i, j])
      end
    end, 1:n-1) |> rm_empty
end


rm_empty(xs::AbstractVector) = map(x -> x, filter(!isnothing, xs))



# function intersect_x(line::Line{T}, x::T) where {T}
#   p0 = line_start(line)
#   x0, y0 = p0.x, p0.y
#   y = k * (x - x0) + y0
#   return Point(x, y)
# end
# 
# function intersect_y(line::Line{T}, y::T) where {T}
#   p0 = line_start(line)
#   x0, y0 = p0.x, p0.y
#   x = (y - y0) / k + x0
#   return Point(x, y)
# end

# function _cellij(ra::SpatRaster, p1::Point, p2::Point)
#   x = (p1.x + p2.x) / 2
#   y = (p1.y + p2.y) / 2

#   b = st_bbox(ra)
#   cellx, celly = st_cellsize(ra)
#   nx, ny = size(ra)[1:2]
#   _location_fast((x, y); b, cellx, celly, nx, ny) # (i, j)
# end
