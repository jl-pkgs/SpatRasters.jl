# function slice2(x::AbstractArray, i, j)
#   cols = repeat([:], ndims(x) - 2)
#   @views x[i, j, cols...]
# end

function resample_first(r::AbstractArray; fact=10, deepcopy=false)
  cols = repeat([:], ndims(r) - 2)

  if deepcopy
    r[1:fact:end, 1:fact:end, cols...]
  else
    @views r[1:fact:end, 1:fact:end, cols...]
  end
end

# for Raster
# st_resample(x::AbstractArray; fact=10, kw...) = resample2(x; fact, kw...)

# function st_resample(ra::SpatRaster; fact=10, deepcopy=false)
#   A = resample2(ra.A; fact, deepcopy)
#   SpatRaster(A, ra.b)
# end

# function st_resample(z::ZArray; fact=10, missingval=0)
#   dat = resample2(z; fact)
#   Raster(dat, st_bbox(z); missingval)
# end

# function st_resample(zs::Vector{<:ZArray}; fact=10, missingval=0)
#   res = Vector{Raster}(undef, length(zs))
#   Threads.@threads for i = eachindex(zs)
#     println("running $i")
#     z = zs[i]
#     dat = st_resample(z; fact)
#     b = st_bbox(z)
#     res[i] = Raster(dat, b)
#   end
#   st_mosaic(res; missingval)
# end

export resample_first
