# using DimensionalData

struct FileNetCDF
  file
end
# struct FileTiff
#   file
# end
struct FileGDAL
  file
end

FileType = Dict(
  ".nc" => FileNetCDF,
  ".tif" => FileGDAL
)

function guess_filetype(f::String)
  ext = file_ext(f)
  FileType[ext](f)
end


st_dims(r::AbstractSpatRaster) = r.lon, r.lat

function st_dims(f::String)
  x = guess_filetype(f)
  st_dims(x)
end


function st_dims(x::FileGDAL)
  gdalinfo(x.file)["dims"]
end

# function st_dims(ra::AbstractRaster)
#   ds = map(dims(ra)) do d
#     DimensionalData.maybeshiftlocus(Center(), d)
#   end
#   x = ds[1].val.data
#   y = ds[2].val.data
#   x, y
# end

st_cellsize(r::AbstractSpatRaster) = st_cellsize(r.lon, r.lat)

function st_cellsize(lon::AbstractVector, lat::AbstractVector)
  median(diff(lon)), median(diff(lat)) # cellx, celly
end

st_cellsize(f::AbstractString) = gdalinfo(f)["cellsize"]

# function st_cellsize(ra::AbstractRaster)
#   # x, y = st_dims(r)
#   # x[2] - x[1], y[2] - y[1]
#   lon = ra.dims[1] # X
#   lat = ra.dims[2] # Y
#   lon[2] - lon[1], lat[2] - lat[1] # cellx, celly
# end
