"""
    st_mosaic(rs::Vector{SpatRaster{T,N}}; kw...) where {T,N}

# Arguments
- `kw`: others to `SpatRaster`
"""
function st_mosaic(rs::Vector{SpatRaster{T,N}}; kw...) where {T,N}
  ra = rs[1]
  # missingval = T(missingval)
  cellsize = st_cellsize(ra)
  box = st_bbox(st_bbox.(rs))
  lon2, lat2 = bbox2dims(box; cellsize)

  _size = length(lon2), length(lat2)
  nd = ndims(ra)
  cols = repeat([:], nd - 2)

  nd >= 3 && (_size = (_size..., size(ra)[3:end]...))
  A = zeros(T, _size)

  # 这里不能使用并行
  # Threads.@threads 
  for i in eachindex(rs)
    ra = rs[i]
    b = st_bbox(ra)
    ilon, ilat = bbox_overlap(b, box; cellsize)
    A[ilon, ilat, cols...] .= ra.A
  end

  (; bands, time, name) = ra
  SpatRaster(A, box; bands, time, name, kw...)
end


function merge_var!(R::AbstractArray, f;
  var=nothing, box::bbox, cellsize=nothing)
  isnothing(cellsize) && (cellsize = st_cellsize(f))
  
  b = st_bbox(f)
  bands = bandnames(f)
  nband = length(bands)
  inds_var = !isnothing(var) ? grep(bands, var) : 1:nband
  ntime = length(inds_var)

  println("Reading data ...")
  A = read_gdal(f, inds_var)
  ilon, ilat = bbox_overlap(b, box; cellsize)
  R[ilon, ilat, 1:ntime] .= A
  nothing
end

"""
    merge_var(fs; vars=nothing, var=nothing,
        progress=true, box::Union{bbox,Nothing}=nothing)

# Arguments
- `fs`: tiff file paths
- `vars`: variables in the tiff
- `var`: variable to merge

# Examples
```julia
box = bbox(-180, -60, 180, 90)
```
"""
function merge_var(fs; vars=nothing, var=nothing,
  progress=true, box::Union{bbox,Nothing}=nothing)
  isnothing(box) && (box = st_bbox(st_bbox.(fs)))

  f = fs[1]
  cellsize = st_cellsize(f)
  lon, lat = bbox2dims(box; cellsize)
  nlon, nlat = length(lon), length(lat)

  bands = bandnames(f)
  ntime = isnothing(vars) ? length(bands) : length(bands) / length(vars)

  R = zeros(Float32, nlon, nlat, ntime)
  p = Progress(ntime)
  for f in fs
    progress && next!(p)
    merge_var!(R, f; var, box, cellsize)
  end
  R
end

export merge_var, merge_var!
