"""
    earth_dist(x1::Matrix, x2::Matrix; R=6378.388)
    earth_dist(p1::Tuple{FT,FT}, x2::AbstractMatrix; R=6378.388) where {FT}

```julia
p1 = [110 30; 111 31]
p2 = [113 32; 115 35]
earth_dist(p1, p2)
```
"""
function earth_dist(x1::Matrix, x2::AbstractMatrix; R=6378.388)
  lon1 = deg2rad.(x1[:, 1])
  lat1 = deg2rad.(x1[:, 2])
  lon2 = deg2rad.(x2[:, 1])
  lat2 = deg2rad.(x2[:, 2])

  _x = @. [cos(lat1) * cos(lon1) cos(lat1) * sin(lon1) sin(lat1)]
  _y = @. [cos(lat2) * cos(lon2) cos(lat2) * sin(lon2) sin(lat2)]
  pp = _x * _y'
  return R .* acos.(clamp.(pp, -1, 1))
end

# 为了追求最大性能, 角度提前转为radian, [in radian]
function earth_dist(lon1::T, lat1::T, lon2::T, lat2::T; R::T=T(6378.388), in_radian=false) where {T}
  if !in_radian
    lon1 = deg2rad(lon1)
    lat1 = deg2rad(lat1)
    lon2 = deg2rad(lon2)
    lat2 = deg2rad(lat2)
  end
  pp = cos(lat1) * cos(lon1) * cos(lat2) * cos(lon2) +
       cos(lat1) * sin(lon1) * cos(lat2) * sin(lon2) +
       sin(lat1) * sin(lat2)
  return R * acos(clamp(pp, -1, 1))
end

# 1 -> 1
function earth_dist(p1::Tuple{T,T}, p2::Tuple{T,T}; R::T=T(6378.388), in_radian=false) where {T<:Real}
  earth_dist(p1[1], p1[2], p2[1], p2[2]; R, in_radian)
end

# 1 -> 1
function earth_dist(p1::AbstractPoint{T}, p2::AbstractPoint{T}; in_radian=false) where {T}
  earth_dist(p1.x, p1.y, p2.x, p2.y; in_radian)
end

# 1 -> n
function earth_dist(p1::Tuple{FT,FT}, x2::AbstractMatrix; R=6378.388) where {FT}
  x1 = [p1[1] p1[2]]
  earth_dist(x1, x2; R)[:]
end


# # 单点
# function distance_earth(p1::AbstractArray{FT}, p2::AbstractArray{FT}; R=6378.388) where {FT}
#   lon1 = deg2rad(p1[1])
#   lat1 = deg2rad(p1[2])
#   lon2 = deg2rad(p2[1])
#   lat2 = deg2rad(p2[2])
#   pp = cos(lat1) * cos(lon1) * cos(lat2) * cos(lon2) +
#        cos(lat1) * sin(lon1) * cos(lat2) * sin(lon2) +
#        sin(lat1) * sin(lat2)
#   return R * acos(clamp(pp, -1, 1))
# end

export earth_dist
