# file = "clusterIds_temporal-(perc_50%,1980-2015).tif"
"""
  $(TYPEDSIGNATURES)

Creates vector polygons for all connected regions of pixels in the raster
sharing a common pixel value.

# Usage

$(METHODLIST)


# Arguments
- `options`: other parameters to 

- `diag` (default `false`): `false` (4 connected) or `true` (8 connected).

- `bands` (default nothing): , bands to be proceeded (e.g, 1:n). If `nothing`,
  it is all bands.

- `drive`: "ESRI Shapefile", "OpenFileGDB", "GPKG"

- `mask`: `img_mask = mask ? band : C_NULL`. All pixels in the mask band with a
  value other than zero will be considered suitable for collection as polygons.


# References

- https://gdal.org/programs/gdal_polygonize.html

- https://github.com/JuliaGeo/GDAL.jl/blob/master/test/drivers.jl

- https://gdal.org/api/gdal_alg.html#_CPPv414GDALPolygonize15GDALRasterBandH15GDALRasterBandH9OGRLayerHiPPc16GDALProgressFuncPv
"""
function gdal_polygonize(raster_file, fout="out.gdb";
  options=String[], diag=false,
  bands=nothing,
  fieldname="grid", nodata=NaN, mask=false, verbose=false, overwrite=true)

  if overwrite && (isfile(fout) || isdir(fout))
    rm(fout; recursive=true)
  end
  if diag
    options = [options..., "8CONNECTED=8"]
  end

  ds_raster = GDAL.gdalopen(raster_file, GDAL.GA_ReadOnly)
  n = nband(raster_file)
  bands === nothing && (bands = 1:n)

  # drive = GDAL.gdalgetdriverbyname("ESRI Shapefile")
  driver = find_shortname(fout)
  _drive = GDAL.gdalgetdriverbyname(driver)
  ds_shp = GDAL.gdalcreate(_drive, fout, 0, 0, 0, GDAL.GDT_Unknown, C_NULL)

  # REF = "WGS_1984"
  REF = GDAL.osrnewspatialreference(C_NULL)
  GDAL.osrimportfromepsg(REF, 4326)

  for i = bands
    verbose && println(i)
    layername = "b_$i"
    band = GDAL.gdalgetrasterband(ds_raster, i)

    if !isnan(nodata)
      GDAL.gdalsetrasternodatavalue(band, nodata)
    end
    # dst_ds = GDAL.ogr_dr_createdatasource(drv, dst_filename, C_NULL)  
    layer = GDAL.gdaldatasetcreatelayer(ds_shp, layername, REF, GDAL.wkbPolygon, C_NULL)

    fielddefn = GDAL.ogr_fld_create(fieldname, GDAL.OFTInteger)
    field = GDAL.ogr_l_createfield(layer, fielddefn, GDAL.TRUE)

    img_mask = mask ? band : C_NULL

    GDAL.gdalpolygonize(
      band,
      img_mask,
      layer, field,
      options, C_NULL, C_NULL)
  end

  GDAL.gdalclose(ds_shp)
  GDAL.gdalclose(ds_raster)
  nothing
end
