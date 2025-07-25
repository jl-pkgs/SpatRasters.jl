function earth_dist(p1::Point3{T}, p2::Point3{T}) where {T}
  earth_dist((p1.x, p1.y), (p2.x, p2.y))
end

"slope in radian"
function dem_slope(p0::Point3{T}, p1::Point3{T}) where {T}
  dl = earth_dist(p0, p1) * 1000 # 水平面上的距离, [km] to [m]
  dz = p1.z - p0.z # [m]
  atan(dz / dl) # radians
end

function dem_slope(p0::Point3{T}, Points::Vector{Point3{T}}) where {T}
  map(p1 -> dem_slope(p0, p1), Points) # αs, H = pi/2 - maximum(αs)
end
