export _st_location, st_location, st_location_exact


"""
    st_location(r::Raster, points::Vector{Tuple{T,T}})

return the overlaping indexes `inds`, and corresponding (i,j)

# Examples
```julia
inds, locs = st_location(r, points)
```
"""
function _st_location(rastsize::RasterSize, x::T, y::T) where {T<:Real}
  (; b, cellsize, nx, ny) = rastsize
  cellx, celly = cellsize
  i = (x - b.xmin) / cellx
  if celly > 0
    j = (y - b.ymin) / celly
  else
    j = (b.ymax - y) / -celly
  end
  i2 = floor(Int, i) + 1
  j2 = ceil(Int, j)
  # @show x, y, i, j, i2, j2 # 存在的问题就是最后一个匹配不上
  if (i2 < 1 || i2 > nx) || (j2 < 1 || j2 > ny)
    return nothing ## TODO, unify return type
  else
    i2, j2
  end
end

function st_location(rastersize::RasterSize, points::Vector{P}; rm_empty::Bool=false) where {
  T<:Real,P<:Union{Tuple{T,T},AbstractPoint{T}}}
  map(p -> _st_location(rastersize, get_x(p), get_y(p)), points)
  # return rm_empty ? _rm_empty(locs) : locs
end


## 全部适用于raster
function st_location(ra::AbstractSpatRaster, args...; kw...)
  rastersize = RasterSize(ra)
  st_location(rastersize, args...; kw...)
end


# - `tol`: (0.5 + tol) * cellsize
function st_location_exact(lon::AbstractVector, lat::AbstractVector, points::Vector{Tuple{T,T}};
  rm_empty::Bool=false, cellsize=nothing, tol=1e-2) where {T<:Real}
  isnothing(cellsize) && (cellsize = st_cellsize(lon, lat))
  cellx, celly = cellsize
  
  map(p -> findnear(p, lon, lat; cellx, celly, tol), points)
  # return rm_empty ? _rm_empty(locs) : locs
end

# function xy2ij(x::T, y::T, b::bbox, cellsize) where {T<:Real}
#   cellx, celly = cellsize
#   i = floor(Int, (x - b.xmin) / cellx)
#   if celly < 0
#     j = floor(Int, (b.ymax - y) / abs(celly))
#   else
#     j = floor(Int, (y - b.ymin) / celly)
#   end
#   return i, j
# end

# function xy2ij(point::Tuple{T,T}, ra::SpatRaster) where {T<:Real}
#   x, y = point
#   xy2ij(x, y, st_bbox(ra), st_cellsize(ra))
# end
