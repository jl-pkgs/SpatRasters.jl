module SpatRasters

using DocStringExtensions
using ProgressMeter
using Statistics: median, cor
using Parameters

export cor, median

## add test data
const dir_data = "$(@__DIR__)/../data" |> abspath
const guanshan_dem = "$dir_data/GuanShan_dem250m.tif"
const guanshan_flowdir_cpp = "$dir_data/GuanShan_flowdir_cpp.tif"
const guanshan_flowdir_gis = "$dir_data/GuanShan_flowdir_gis.tif"
export guanshan_dem, guanshan_flowdir_cpp, guanshan_flowdir_gis

export bbox, in_bbox, bbox_overlap
export bbox2lims,
  bbox2cellsize,
  bbox2range, bbox2vec,
  bbox2dims, bbox2ndim
export range2bbox
export st_points
export st_bbox, st_dims, st_cellsize, st_mosaic
export st_write, st_read, nlyr
export rm_shp
export getgeotransform
export read_sf, write_sf

export xy2ij, cellArea

export write_gdal, read_gdal
export ogr_info, gdal_info, gdalinfo, gdal_nodata
export bandnames, set_bandnames, nband, nlayer
export gdal_polygonize

function read_sf end
function write_sf end

function nband end
function nlayer end
function gdal_polygonize end
function read_gdal end
function write_gdal end
function gdalinfo end
function gdal_info end
function ogr_info end
function bandnames end
function set_bandnames end
function gdal_nodata end

nlyr = nband
st_write = write_gdal
st_read = read_gdal

include("tools_Ipaper.jl")
include("bbox.jl")
include("SpatRaster.jl")
include("Ops.jl")

include("st_bbox.jl")
include("st_dims.jl")
include("IO.jl")

include("st_extract.jl")
include("st_resample.jl")
include("st_mosaic.jl")
include("st_crop.jl")

include("distance.jl")
include("hydro/Hydro.jl")

include("Interp/Interp.jl")
include("methods/intersect.jl")
include("terrain/terrain.jl")

export st_coords

function st_coords(ra::SpatRaster)
  lon, lat = st_dims(ra)
  Lon, Lat = meshgrid(lon, lat)
  [Lon[:] Lat[:]]
end

function st_points(x::AbstractVector, y::AbstractVector)
  [(x[i], y[i]) for i in eachindex(x)]
end


function shp_files(f)
  [f,
    replace(f, ".shp" => ".shx"),
    replace(f, ".shp" => ".prj"),
    replace(f, ".shp" => ".dbf")]
end

function rm_shp(f)
  rm.(shp_files(f))
  nothing
end


function xy2ij(x::T, y::T, b::bbox, cellsize) where {T<:Real}
  cellx, celly = cellsize
  i = floor(Int, (x - b.xmin) / cellx)
  if celly < 0
    j = floor(Int, (b.ymax - y) / abs(celly))
  else
    j = floor(Int, (y - b.ymin) / celly)
  end
  return i, j
end

function xy2ij(point::Tuple{T,T}, ra::SpatRaster) where {T<:Real}
  x, y = point
  xy2ij(x, y, st_bbox(ra), st_cellsize(ra))
end

# in km^2
function cellArea(x, y, cellsize)
  cellx, celly = cellsize
  dx = earth_dist((x, y), (x + cellx, y))
  dy = earth_dist((x, y), (x, y + celly))
  return dx * dy
end

end # module SpatRasters
