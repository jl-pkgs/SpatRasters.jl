export weight_idw!, weight_adw!


function weight_idw!(neighbor::Neighbor{FT,N}; m::Int=2) where {FT,N}
  (; count, dist, weight, dims) = neighbor
  weight .= FT(0.0)
  nlon, nlat = dims[1:2]
  for i in 1:nlon, j in 1:nlat
    n = count[i, j]
    # _w = @view weight[i, j, 1:n]
    _dist = @view dist[i, j, 1:n]
    wk = @.(1 / _dist^m)
    wk .= wk ./ sum(wk)
    weight[i, j, 1:n] .= wk
  end
end


function weight_adw!(neighbor::Neighbor{FT,N};
  cdd=100, m::Int=2) where {FT,N}

  (; nmax, count, dims, weight, dist) = neighbor
  angle::Array{FT,N} = neighbor.angle

  weight .= FT(0.0)
  nlon, nlat = dims[1:2]
  ∅ = FT(0)
  wks = [zeros(nmax) for _ in 1:Threads.nthreads()]
  αs = [zeros(nmax) for _ in 1:Threads.nthreads()]

  @inbounds @threads for j in 1:nlon
    for i in 1:nlat
      wk = wks[Threads.threadid()]
      α = αs[Threads.threadid()]
      @. wk = ∅
      @. α = ∅
      n = count[i, j]
      n == 0 && continue

      for k in 1:n
        wk[k] = exp(-dist[i, j, k] / cdd)^m
      end

      for k in 1:n # for each candidates, update `w` according to `angle`
        ∑ = ∑w = ∅
        @inbounds for l in 1:n
          k == l && continue
          Δθ::FT = deg2rad(angle[i, j, k] - angle[i, j, l])
          w1 = wk[l]
          ∑ += w1 * FT(1 - cos(Δθ))  # Xavier 2016, Eq. 7
          ∑w += w1                   # Xavier 2016, Eq. 8
        end
        α[k] = ∑ / ∑w
      end # end of candidates

      @. wk *= (1 + α)               # Xavier 2016, Eq. 9
      ∑w = sum(wk)
      @. wk /= ∑w
      @. weight[i, j, 1:n] = wk[1:n]
    end # end of lon
  end # end of lat
end

# used for test
function _weight_adw(dist::AbstractVector{FT}, angle::AbstractVector{FT}; cdd::FT=FT(100), m::Int=2) where {FT}
  n = length(dist)
  wk = zeros(FT, n)
  α = zeros(FT, n)
  for k in 1:n
    wk[k] = exp(-dist[k] / cdd)^m
  end
  ∅ = FT(0)
  for k in 1:n # for each candidates, update `w` according to `angle`
    ∑ = ∑w = ∅
    @inbounds for l in 1:n
      k == l && continue
      Δθ::FT = deg2rad(angle[k] - angle[l])
      w1 = wk[l]
      ∑ += w1 * FT(1 - cos(Δθ))  # Xavier 2016, Eq. 7
      ∑w += w1                   # Xavier 2016, Eq. 8
    end
    α[k] = ∑ / ∑w
  end # end of candidates
  @. wk *= (1 + α)               # Xavier 2016, Eq. 9
  ∑w = sum(wk)
  @. wk /= ∑w
  wk
end

export _weight_adw
