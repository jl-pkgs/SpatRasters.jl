"""
    bbox(xmin, ymin, xmax, ymax)
    bbox(;xmin, ymin, xmax, ymax)
    bbox2tuple(b::bbox)
    bbox2vec(b::bbox)
    bbox2affine(size::Tuple{Integer, Integer}, b::bbox)

Spatial bounding box
"""
Base.@kwdef struct bbox
  xmin::Float64
  ymin::Float64
  xmax::Float64
  ymax::Float64
end


# size: array size
function bbox2cellsize(b::bbox, size)
  nlon, nlat = size[1:2]
  cellx = (b.xmax - b.xmin) / nlon
  celly = (b.ymax - b.ymin) / nlat
  cellx, -celly # 默认y是倒序
end

# 两个至少提供一个
"""
    bbox2dims(b::bbox; size=nothing, cellsize=nothing, reverse_lat=true)

# Arguments
- `cellsize`: celly为负，则lat倒序
"""
function bbox2dims(b::bbox; size=nothing, cellsize=nothing, reverse_lat=true)
  if size !== nothing && cellsize === nothing
    cellsize = bbox2cellsize(b, size)
  end
  length(cellsize) == 1 && (cellsize = [1, 1] .* cellsize)

  cellx, celly = abs.(cellsize)
  lon = b.xmin+cellx/2:cellx:b.xmax
  lat = b.ymin+celly/2:celly:b.ymax

  (cellsize[2] < 0 || reverse_lat) && (lat = reverse(lat))
  lon, lat
end


function bbox2ndim(b; size=nothing, cellsize=nothing)
  x, y = bbox2dims(b; size, cellsize)
  length(x), length(y)
end

"""
    in_bbox(b::bbox, (lon, lat))  
    in_bbox(bs::Vector{bbox}, (lon, lat))
"""
in_bbox(b::bbox, (lon, lat)) = (b.xmin < lon < b.xmax) && (b.ymin < lat < b.ymax)

in_bbox(bs::Vector{bbox}, (lon, lat)) = [in_bbox(b, (lon, lat)) for b in bs]


bbox2range(b::bbox) = [b.xmin, b.xmax, b.ymin, b.ymax]
bbox2tuple(b::bbox) = (xmin=b.xmin, ymin=b.ymin, xmax=b.xmax, ymax=b.ymax)
bbox2vec(b::bbox) = [b.xmin, b.ymin, b.xmax, b.ymax]
bbox2lims(b::bbox) = ((b.xmin, b.xmax), (b.ymin, b.ymax))

range2bbox(r::AbstractVector) = bbox(r[1], r[3], r[2], r[4])

function bbox_overlap(b::bbox, box::bbox; size=nothing, cellsize=nothing, reverse_lat=true, zip=true)
  Lon, Lat = bbox2dims(box; size, cellsize, reverse_lat)

  ilon = findall(b.xmin .< Lon .< b.xmax) |> _zip
  ilat = findall(b.ymin .< Lat .< b.ymax) |> _zip

  if zip
    ilon = _zip(ilon)
    ilat = _zip(ilat)
  end

  ## 若`b`的范围 > `box`，则会出现不匹配的现象，因此下述检查不能开启
  # lon, lat = bbox2dims(b; size, cellsize, reverse_lat)
  # @assert length(lon) == length(ilon)
  # @assert length(lat) == length(ilat)
  ilon, ilat
end

_zip(x) = x[1]:x[end]
