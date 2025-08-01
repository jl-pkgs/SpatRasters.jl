export find_neighbor, Neighbor
export distance_norm

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


include("angle.jl")
include("find_neighbor.jl")
include("weights.jl")
include("interp_tps.jl")
include("bilinear_helper.jl")
include("bilinear.jl")
include("bilinear/bilinear_irregular.jl")


"""
    interp(ra::SpatRaster, target::SpatRaster;
        nmax::Int=20, radius::Real=100, do_angle=false, method="adw", kw...)
    
    interp(x::AbstractMatrix, y::AbstractArray{FT}, target::SpatRaster;
        nmax::Int=20, radius::Real=200, do_angle=false, method="adw", kw...)
    
    interp(neighbor::Neighbor{FT,N}, y::AbstractArray{FT}, target::SpatRaster; ignored...)

## Arguments
- `method`: one of "adw", "idw", "nearest", "bilinear"

- `kw`: other parameters to `wfun`
  + `weight_idw!`: `m`
  + `weight_adw!`: `cdd`, `m`
"""
function interp(ra::SpatRaster, target::SpatRaster;
  nmax::Int=20, radius::Real=100, do_angle=false, method="adw", kw...)
  dims = size(ra)
  nlon, nlat = dims[1:2]
  ntime = length(dims) >= 3 ? dims[3] : 1

  if method == "bilinear"
    x, y = st_dims(ra)
    xx, yy = st_dims(target)
    Z = ra.A
    R = bilinear(x, y, Z, xx, yy)
    rast(R, st_bbox(target))
  else
    X = st_coords(ra) # [lon, lat]
    Y = reshape(ra.A, nlon * nlat, ntime)

    interp(X, Y, target; nmax, radius, do_angle, method, kw...)
  end
end


function interp(x::AbstractMatrix, y::AbstractArray{FT}, target::SpatRaster;
  nmax::Int=20, radius::Real=100, do_angle=false,
  method="adw", kw...) where {FT}

  neighbor = find_neighbor(target, x; nmax, radius, do_angle)

  if method == "nearest"
    interp_nearest(neighbor, y, target)
  else
    wfun = WFUNS[method]
    wfun(neighbor; kw...)
    interp_weight(neighbor, y, target)
  end
end


WFUNS = Dict(
  "idw" => weight_idw!,
  "adw" => weight_adw!,
)


function interp_weight(neighbor::Neighbor{FT,N}, y::AbstractArray{FT}, target::SpatRaster; ignored...) where {FT,N}
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


function interp_nearest(neighbor::Neighbor{FT,N}, y::AbstractArray{FT}, target::SpatRaster; ignored...) where {FT,N}
  ntime = size(y, 2)
  lon, lat = st_dims(target)
  nlon, nlat = length(lon), length(lat)
  R = zeros(FT, nlon, nlat, ntime)

  (; count, index, weight) = neighbor
  # ∅ = FT(0)
  @inbounds for i in 1:nlon
    for j in 1:nlat
      n_control = count[i, j]
      inds = @view index[i, j, 1:n_control]

      for k in 1:ntime
        for l in 1:n_control # control points
          _i = inds[l]
          yᵢ = y[_i, k]
          notnan = yᵢ == yᵢ
          if notnan # 最邻近原则，只要不是Nan，就可以跳出
            R[i, j, k] = yᵢ
            break # 跳出for循环
          end
        end # end n_control
      end # end time
    end # end lat
  end # end lon
  rast(R, target)
end


export interp
