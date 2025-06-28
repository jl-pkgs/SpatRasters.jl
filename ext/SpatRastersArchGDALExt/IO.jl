function read_gdal(file::AbstractString, options...)
  ArchGDAL.read(file) do dataset
    ArchGDAL.read(dataset, options...)
  end
end

# convert into a rast
function read_gdal(file::AbstractString, indices, b::bbox, options...)
  # cellsize = 1 / 1200
  cellsize = gdalinfo(file)["cellsize"]
  box = st_bbox(file)
  lon, lat = st_dims(file)
  ilon, ilat = bbox_overlap(b, box; cellsize)

  _lon, _lat = lon[ilon], lat[ilat]
  _b = st_bbox(_lon, _lat)
  _A = read_gdal(file, indices, ilat, ilon, options...)
  rast(_A, _b; nodata=gdal_nodata(file))
end

function read_gdal(file::AbstractString, b::bbox, options...)
  indices = 1:nband(file)
  read_gdal(file, indices, b, options...)
end

## This part is borrowed from the GeoArrays.jl package.
# MIT License, Copyright (c) 2018 Maarten Pronk
# <https://github.com/evetion/GeoArrays.jl/blob/master/src/io.jl>

# using GeoFormatTypes, ArchGDAL
# WGS84 = convert(WellKnownText, EPSG(4326))
# GFT.val(ga.crs)
const WGS84 = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]"

# const OPTIONS_DEFAULT_TIFF = Dict(
#   "TILED" => "YES", # not work
#   "COMPRESS" => "DEFLATE"
# )
import ArchGDAL: OF_UPDATE

function gdal_setproj!(f::AbstractString, transform::Vector{Cdouble})
  # ArchGDAL.open(f, "r+") do ds
  ArchGDAL.read(f; flags=OF_UPDATE) do ds
    ## Set geotransform and crs
    ArchGDAL.GDAL.gdalsetgeotransform(ds.ptr, transform)
    ArchGDAL.GDAL.gdalsetprojection(ds.ptr, WGS84)
  end
  return nothing
end


function write_gdal(data::AbstractArray, f::AbstractString;
  nodata=nothing, options=String[], NUM_THREADS=4, BIGTIFF=false)

  dtype = eltype(data)
  shortname = find_shortname(f)
  driver = ArchGDAL.getdriver(shortname)

  width, height = size(data)[1:2]
  ndims(data) == 2 && (nbands = 1)
  ndims(data) == 3 && (nbands = size(data, 3))

  if !isnothing(nodata) && !isa(nodata, Vector)
    nodata = fill(nodata, nbands)
  end

  if (shortname == "GTiff")
    options = [options..., "COMPRESS=DEFLATE", "TILED=YES", "NUM_THREADS=$NUM_THREADS"]
    BIGTIFF && (push!(options, "BIGTIFF=YES"))
  end

  try
    convert(ArchGDAL.GDALDataType, dtype)
  catch
    dtype, data = cast_to_gdal(data)
  end

  ArchGDAL.create(f; driver, width, height, nbands, dtype, options) do dataset
    for i = 1:nbands
      band = ArchGDAL.getband(dataset, i)
      ArchGDAL.write!(band, data[:, :, i])
      !isnothing(nodata) && ArchGDAL.GDAL.gdalsetrasternodatavalue(band.ptr, nodata[i])
    end
  end
end

# only support WGS84 proj
function write_gdal(ra::AbstractSpatRaster, f::AbstractString;
  nodata=nothing, options=String[], NUM_THREADS=4, BIGTIFF=true)

  isnothing(nodata) && (nodata = ra.nodata)
  write_gdal(ra.A, f; nodata, options, NUM_THREADS, BIGTIFF)
  gdal_setproj!(f, getgeotransform(ra))

  !isnothing(ra.bands) && set_bandnames(f, ra.bands)
  return f
end

# # Slice data and replace missing by nodata
# if isa(dtype, Union) && dtype.a == Missing
#   dtype = dtype.b
#   try
#     convert(ArchGDAL.GDALDataType, dtype)
#   catch
#     dtype, data = cast_to_gdal(data)
#   end
#   nodata === nothing && (nodata = typemax(dtype))
#   m = ismissing.(data)
#   data[m] .= nodata
#   data = Array{dtype}(data)
#   use_nodata = true
# end
