# using ArchGDAL
# GDAL.gdalgetgeotransform(ds)
function getgeotransform!(
  dataset::Ptr,
  transform::Vector{Cdouble},
)::Vector{Cdouble}
  @assert length(transform) == 6
  result = GDAL.gdalgetgeotransform(dataset, pointer(transform))
  if result != GDAL.CE_None
    # The default geotransform.
    transform .= (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
  end
  return transform
end

getgeotransform(dataset::Ptr)::Vector{Cdouble} =
  getgeotransform!(dataset, Vector{Cdouble}(undef, 6))

width(band) = GDAL.gdalgetrasterbandxsize(band) # error at here
height(band) = GDAL.gdalgetrasterbandysize(band)

function gdalinfo(file::AbstractString)
  # nband = nband(file)
  ds = gdal_open(file)
  gt = getgeotransform(ds)
  band = GDAL.gdalgetrasterband(ds, 1)
  w, h = width(band), height(band)
  gdal_close(ds)

  dx, dy = gt[2], gt[end]
  x0 = gt[1] #+ dx/2
  y0 = gt[4] #- dy/2
  x1 = x0 + w * dx
  y1 = y0 + h * dy

  lon = x0+dx/2:dx:x1
  lat = y0+dy/2:dy:y1
  b = st_bbox(lon, lat) # b = bbox(x0, y0, x1, y1)

  Dict(
    "file" => basename(file),
    "bbox" => b,
    "cellsize" => [dx, dy],
    "dims" => [lon, lat],
    "size" => Int64.([w, h, nband(file)])
  ) # convert to tuple
end
