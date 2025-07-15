export weight_idw!

@with_kw struct Neighbor{FT,N}
  nmax::Int = 20
  dims::Tuple = ()
  count::AbstractArray{Integer} = zeros(Int, dims...)
  index::AbstractArray{Integer,N} = zeros(Int, dims..., nmax)
  dist::AbstractArray{FT,N} = zeros(FT, dims..., nmax)
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

function find_neighbor(ra::SpatRaster, points_source; nmax::Int=20, radius::Real=200)
  X = collect(points_source')
  tree = BallTree(X, Haversine(6371.0); reorder=false)

  lon, lat = st_dims(ra)
  nlon, nlat = length(lon), length(lat)
  neighbor = Neighbor(nmax, (nlon, nlat))
  (; count, index, dist) = neighbor

  for i in 1:nlon, j in 1:nlat
    p_target = [lon[i], lat[j]]

    _inds, _dists = knn(tree, p_target, nmax)
    _I = _dists .<= radius
    _inds = @view _inds[_I]
    _dists = @view _dists[_I]

    n = length(_dists)
    if n > 0
      count[i, j] = n
      index[i, j, 1:n] .= _inds
      dist[i, j, 1:n] .= _dists
    end
  end
  neighbor
end


function weight_idw!(neighbor::Neighbor{FT,N}; m=2) where {FT,N}
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

