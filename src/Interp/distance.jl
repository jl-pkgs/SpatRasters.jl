function distance_norm(x::P, y::P) where {P}
  @assert length(x) == length(y)
  dist2 = 0.0
  @inbounds for i in eachindex(x)
    dist2 += (x[i] - y[i])^2
  end
  sqrt(dist2)
end

