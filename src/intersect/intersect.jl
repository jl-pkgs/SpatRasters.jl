include("datatype.jl")
import Base: intersect

export intersect
export line_start, line_end
export azimuth2slope
export SVF_azimuth, SVF_azimuth_simple
export SVF


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


"和x轴的交点，如果k = Inf or -Inf，则需要用y轴相交的方法"
function intersect_x(line::Line{T}, xs::AbstractVector{T}) where {T}
  (; k) = line
  p0 = line_start(line)
  x0, y0 = p0.x, p0.y
  map(x -> begin
      y = k * (x - x0) + y0
      Point(x, y)
    end, xs)
end

function intersect_y(line::Line{T}, ys::AbstractVector{T}) where {T}
  (; k) = line
  p0 = line_start(line)
  x0, y0 = p0.x, p0.y
  map(y -> begin
      x = (y - y0) / k + x0
      Point(x, y)
    end, ys)
end


function range_lat(values::AbstractVector{T}, min::T, max::T) where {T}
  ibeg = searchsortedfirst(values, max) # 逆序
  iend = searchsortedlast(values, min)
  ibeg:iend
end

function range_lon(values::AbstractVector{T}, min::T, max::T) where {T}
  ibeg = searchsortedfirst(values, min)
  iend = searchsortedlast(values, max)
  ibeg:iend
end


function intersect(ra::SpatRaster, line::Line; cellsize=nothing)
  isnothing(cellsize) && (cellsize = st_cellsize(ra))
  cellx, celly = cellsize

  b = st_bbox(ra)
  lon = b.xmin:cellx:b.xmax # 采用的是网格边界
  lat = celly < 0 ? (b.ymax:celly:b.ymin) : (b.ymin:celly:b.ymax)

  bl = st_bbox(line)
  ## 判断交点
  ilat = range_lat(lat, bl.ymin, bl.ymax) # bl.ymin .<= lat .<= bl.ymax # 
  ilon = range_lon(lon, bl.xmin, bl.xmax) # bl.xmin .<= lon .<= bl.xmax # 判断大致的范围
  xs = @view lon[ilon]
  ys = @view lat[ilat]

  points_y = intersect_y(line, ys) # 与水平线的交点
  points_x = intersect_x(line, xs) # 与垂线的交点

  points = cat(points_y, points_x, dims=1) |> rm_empty

  x = map(p -> p.x, points)
  inds = sortperm(x)
  points = @view points[inds] # 对points进行排序  

  ## 然后两点判断一个网格位置
  _cellij(ra, points; cellsize)
end


function earth_dist(p1::Point3{T}, p2::Point3{T}) where {T}
  earth_dist((p1.x, p1.y), (p2.x, p2.y))
end


"""
仅用于计算相交

# Arguments
points: 与网格边界相交的所有点

# Return
每两个点 确定一个网格中心，返回的是网格中心的[x, y, elev]
"""
function _cellij(ra::SpatRaster, points::AbstractVector{Point{T}}; cellsize=nothing) where {T}
  isnothing(cellsize) && (cellsize = st_cellsize(ra))

  n = length(points)
  b = st_bbox(ra)
  lon, lat = st_dims(ra)
  cellx, celly = cellsize
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
