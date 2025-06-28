# include("bbox.jl")
st_bbox(r::AbstractSpatRaster) = r.b

function st_bbox(v::Vector{<:Real})
  bbox(v[1:4]...)
end

function st_bbox(lon, lat)
  cellx = abs(lon[2] - lon[1])
  celly = abs(lat[2] - lat[1])
  bbox(minimum(lon) - cellx / 2, minimum(lat) - celly / 2,
    maximum(lon) + cellx / 2, maximum(lat) + celly / 2)
end

function st_bbox(f::String)
  lon, lat = st_dims(f)
  st_bbox(lon, lat)
end


# get large bbox, also know as bbox_merge
function st_bbox(bs::Vector{bbox})
  bboxs = bbox2vec.(bs)
  bbox_mat = cat(bboxs..., dims=2)
  b = [minimum(bbox_mat[1:2, :], dims=2)..., maximum(bbox_mat[3:4, :], dims=2)...]
  bbox(b...)
end

## add support for Rasters
# using Extents
# function st_bbox(ra::AbstractRaster)
#   # x, y = st_dims(ra)
#   xlim, ylim = Extents.extent(ra)
#   bbox(xlim[1], ylim[1], xlim[2], ylim[2])
# end

# st_bbox(ras::Vector{<:Raster}) = st_bbox(st_bbox.(ras))

# st_bbox(z::ZArray) = Terra.bbox(z.attrs["bbox"]...)
# st_bbox(zs::Vector{<:ZArray}) = st_bbox(st_bbox.(zs))


# st_bbox(f::String) = gdalinfo(f)["bbox"]
