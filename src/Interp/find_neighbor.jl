export weight_idw!

@with_kw struct Neighbor{FT,N}
  nmax::Int = 20
  dims::Tuple = ()
  count::AbstractArray{Integer} = zeros(Int, dims...)
  index::AbstractArray{Integer,N} = zeros(Int, dims..., nmax)
  "distance (km)"
  dist::AbstractArray{FT,N} = zeros(FT, dims..., nmax)
  "azimuth angle (degree, 0-360), N: 0°, E: 90°"
  angle::AbstractArray{FT,N} = zeros(FT, dims..., nmax)
  weight::AbstractArray{FT,N} = zeros(FT, dims..., nmax)
end


function Neighbor(nmax::Int, dims::Tuple=(); FT=Float64)
  N = length(dims) + 1
  Neighbor{FT,N}(; nmax, dims)
end


function find_neighbor(point, tree; nmax::Int=20, radius::Real=200, m=2)
  inds, dists = knn(tree, point, nmax)
  I = dists .<= radius
  inds = @view inds[I]
  dists = @view dists[I]
  (; index=inds, distance=dists)
end

function find_neighbor(ra::SpatRaster, X; nmax::Int=20, radius::Real=200,
  do_angle=true)

  Xt = collect(X')
  tree = BallTree(Xt, Haversine(6371.0); reorder=false)

  lon, lat = st_dims(ra)
  nlon, nlat = length(lon), length(lat)
  neighbor = Neighbor(nmax, (nlon, nlat))
  (; count, index, dist, angle) = neighbor

  for i in 1:nlon, j in 1:nlat
    p0 = [lon[i], lat[j]]

    _inds, _dists = knn(tree, p0, nmax)
    _I = sortperm(_dists)
    _I = _I[_dists[_I].<=radius]

    _inds = @view _inds[_I]
    _dists = @view _dists[_I]

    n_control = length(_dists)
    n_control == 0 && continue

    count[i, j] = n_control
    index[i, j, 1:n_control] .= _inds
    dist[i, j, 1:n_control] .= _dists

    if do_angle
      for c in 1:n_control
        _i = _inds[c]
        p1 = @view X[_i, :]
        angle[i, j, c] = angle_azimuth_sphere(p0, p1)
      end
    end
  end
  neighbor
end


function weight_idw!(neighbor::Neighbor{FT,N}; m::Int=2) where {FT,N}
  (; count, dist, weight, dims) = neighbor
  weight .= FT(0.0)
  nlon, nlat = dims[1:2]
  for i in 1:nlon, j in 1:nlat
    n = count[i, j]
    # _w = @view weight[i, j, 1:n]
    _dist = @view dist[i, j, 1:n]
    _w = @.(1 / _dist^m)
    _w .= _w ./ sum(_w)
    weight[i, j, 1:n] .= _w
  end
end

function weight_adw!(neighbor::Neighbor{FT,N}; m::Int=2) where {FT,N}
  (; count, dist, weight, dims) = neighbor
  weight .= FT(0.0)
  nlon, nlat = dims[1:2]
  for i in 1:nlon, j in 1:nlat
    n = count[i, j]
    # _w = @view weight[i, j, 1:n]
    _dist = @view dist[i, j, 1:n]
    _w = @.(1 / _dist^m)
    _w .= _w ./ sum(_w)
    weight[i, j, 1:n] .= _w
  end
end


"""
"""
function interp(x::AbstractMatrix, y::AbstractArray{FT}, target::SpatRaster;
  nmax::Int=20, radius::Real=200, 
  do_angle=false,
  wfun::Function=weight_idw!, kw...) where {FT}

  neighbor = find_neighbor(target, x; nmax, radius, do_angle)
  wfun(neighbor; kw...)

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
