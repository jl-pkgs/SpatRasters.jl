export find_neighbor, Neighbor
export distance_norm, distance_earth

using Base.Threads
using NearestNeighbors
using Distances

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


include("angle.jl")
include("find_neighbor.jl")
include("weights.jl")
include("interp_tps.jl")


"""
    interp(x::AbstractMatrix, y::AbstractArray{FT}, target::SpatRaster;
        nmax::Int=20, radius::Real=200, do_angle=false,
        wfun::Function=weight_idw!, kw...)
    interp(neighbor::Neighbor{FT,N}, y::AbstractArray{FT}, target::SpatRaster; ignored...)

## Arguments
- `kw`: other parameters to `wfun`
  + `weight_idw!`: `m`
  + `weight_adw!`: `cdd`, `m`
"""
function interp(x::AbstractMatrix, y::AbstractArray{FT}, target::SpatRaster;
  nmax::Int=20, radius::Real=200, do_angle=false,
  wfun::Function=weight_idw!, kw...) where {FT}

  neighbor = find_neighbor(target, x; nmax, radius, do_angle)
  wfun(neighbor; kw...)
  interp(neighbor, y, target)
end

function interp(neighbor::Neighbor{FT,N}, y::AbstractArray{FT}, target::SpatRaster; ignored...) where {FT, N}
  ntime = size(y, 2)
  lon, lat = st_dims(target)
  nlon, nlat = length(lon), length(lat)
  R = zeros(FT, nlon, nlat, ntime)

  (; count, index, weight) = neighbor
  ∅ = FT(0)
  for i in 1:nlon, j in 1:nlat
    n_control = count[i, j]
    inds = @view index[i, j, 1:n_control]
    ws = @view weight[i, j, 1:n_control]

    for k in 1:ntime
      ∑ = FT(0)
      ∑w = 0
      for l in 1:n_control # control points
        _i = inds[l]
        yᵢ = y[_i, k]
        notnan = yᵢ == yᵢ
        ∑ += ifelse(notnan, yᵢ * ws[l], ∅)
        ∑w += ifelse(notnan, ws[l], ∅)
      end
      R[i, j, k] = ∑ / ∑w
    end
  end
  rast(R, target)
end

export interp
