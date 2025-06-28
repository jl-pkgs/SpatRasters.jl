export SpatRastersArchGDALExt
module SpatRastersArchGDALExt

export write_tiff

using DocStringExtensions: TYPEDSIGNATURES, METHODLIST

using ArchGDAL
using ArchGDAL.GDAL
using ArchGDAL.GDAL.GDAL_jll: gdalinfo_path, ogrinfo_path

using SpatRasters
import SpatRasters: write_gdal, read_gdal, gdalinfo, getgeotransform, 
  gdal_polygonize, nband, nlayer, gdal_nodata, 
  gdal_info, ogr_info, bandnames, set_bandnames, 
  find_shortname, cast_to_gdal

# import SpatRasters: WGS84

include("gdal_basic.jl")
include("IO.jl")
include("gdalinfo.jl")
include("gdal_polygonize.jl")

export gdal_polygonize
export nband, nlayer
export gdal_nodata
export bandnames, set_bandnames
export gdal_info, ogr_info
export write_gdal, read_gdal

end
