export st_extract

function getvalue(A::AbstractArray, i::Int, j::Int)
  cols = repeat([:], ndims(ra) - 2)
  A[i, j, cols...]
end

## TODO: `_rm_empty`让情况变得复杂
function _rm_empty(x::Vector)
  inds = findall(!isnothing, x)
  inds, x[inds]
end

# function _is_empty(locs::Vector)
#   isempty = map(p -> (p[1] != -1 && p[2] != -1), locs) # 有一个为-1就是空
#   isempty, locs
# end

function st_extract(ra::AbstractSpatRaster{FT}, point::P;) where {
  FT,T,P<:Union{Tuple{T,T},AbstractPoint{T}}}

  ij = st_location(RasterSize(ra), point)
  isnothing(ij) && return FT(NaN)
  getvalue(ra.A, ij...)
end


function st_extract(ra::AbstractSpatRaster, points::Vector{P}; combine=hcat) where {
  T<:Real,P<:Union{Tuple{T,T},AbstractPoint{T}}}

  _locs = st_location(ra, points)
  inds, locs = _rm_empty(_locs)

  cols = repeat([:], ndims(ra) - 2)
  lst = [ra.A[i, j, cols...] for (i, j) in locs]
  (; index=inds, value=combine(lst...))
end
