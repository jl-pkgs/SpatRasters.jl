export angle_azimuth_sphere, angle_azimuth

"""
球面方位角 (0~360°)

```julia
angle_azimuth_sphere((0, 0), (0, 1)) # 正北，0
angle_azimuth((0, 0), (0, 1))        # 正北，0
```

## References
https://chatgpt.com/c/6876ff2c-834c-8012-98cd-2bc9691d453f
"""
function angle_azimuth_sphere(p1::AbstractVector{FT}, p2::AbstractVector{FT}; to_degree=true) where {FT}
  lon1, lat1 = p1
  lon2, lat2 = p2

  φ1 = deg2rad(lat1)
  φ2 = deg2rad(lat2)
  Δλ = deg2rad(lon2 - lon1)

  x = sin(Δλ) * cos(φ2)
  y = cos(φ1) * sin(φ2) - sin(φ1) * cos(φ2) * cos(Δλ)
  θ = atan(x, y)
  !to_degree && return θ
  
  azimuth = rad2deg(θ)
  azimuth < 0 && (azimuth += 360) # 标准化方位角为0-360度
  return azimuth
end

function _radian2_azimuth(θ)
  azimuth = rad2deg(θ)
  azimuth < 0 && (azimuth += 360) # 标准化方位角为0-360度
  return azimuth
end

function angle_azimuth(p1::AbstractVector{FT}, p2::AbstractVector{FT}; to_degree=true) where {FT}
  lon1, lat1 = p1
  lon2, lat2 = p2
  dx = lon2 - lon1
  dy = lat2 - lat1
  # 需要atan2(Δx, Δy)来计算从北向东的方位角
  # azimuth = -rad2deg(atan(dy, dx)) + 90
  !to_degree && return θ
  
  θ = atan(dx, dy)
  azimuth = rad2deg(θ)
  azimuth < 0 && (azimuth += 360) # 标准化方位角为0-360度
  return azimuth
end


function angle_azimuth_sphere(p1::AbstractVector{FT}, p2::AbstractMatrix{FT}) where {FT}
  map(p -> angle_azimuth_sphere(p1, p), eachrow(p2))
end

function angle_azimuth(p1::AbstractVector{FT}, p2::AbstractMatrix{FT}) where {FT}
  map(p -> angle_azimuth(p1, p), eachrow(p2))
end
