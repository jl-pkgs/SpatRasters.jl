export st_crop

st_crop(f::String, b::bbox) = read_gdal(f, b)

function st_crop(ra::SpatRaster, b::bbox)
  (; name, time, nodata, bands) = ra
  box = st_bbox(ra)
  lon, lat = st_dims(ra)
  cellsize = st_cellsize(ra)
  ix, iy = bbox_overlap(b, box; cellsize, reverse_lat=true)
  _lon, _lat = lon[ix], lat[iy]
  _b = st_bbox(_lon, _lat)

  A = ra.A
  cols = repeat([:], ndims(A) - 2)
  _A = A[ix, iy, cols...]
  rast(_A, _b; name, time, nodata, bands)
end
