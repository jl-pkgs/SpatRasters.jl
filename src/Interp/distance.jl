export distance_norm, distance_earth


function distance_norm(p1::AbstractArray{FT}, p2::AbstractArray{FT}; kw...) where {FT}
  # @assert length(p1) == length(p2)
  dist2 = FT(0.0)
  @inbounds for i in eachindex(p1)
    dist2 += (p1[i] - p2[i])^2
  end
  sqrt(dist2)
end

function distance_earth(p1::AbstractArray{FT}, p2::AbstractArray{FT}; R=6378.388) where {FT}
  lon1 = deg2rad(p1[1])
  lat1 = deg2rad(p1[2])
  lon2 = deg2rad(p2[1])
  lat2 = deg2rad(p2[2])

  pp = cos(lat1) * cos(lon1) * cos(lat2) * cos(lon2) +
       cos(lat1) * sin(lon1) * cos(lat2) * sin(lon2) +
       sin(lat1) * sin(lat2)
  return R * acos(clamp(pp, -1, 1))
end
