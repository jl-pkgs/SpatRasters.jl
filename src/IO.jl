export read_gdal

"""
    read_gdal(file::String, options...)
    read_gdal(files::Array{String,1}, options)

# Arguments:
- `options`: other parameters to `ArchGDAL.read(dataset, options...)`.

# Return
"""
# read multiple tiff files and cbind
function read_gdal(files::Vector{<:AbstractString}, options...)
  # bands = collect(bands)
  # bands = collect(Int32, bands)
  res = map(file -> read_gdal(file, options...), files)
  res
  # vcat(res...)
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
# const drivers = AG.listdrivers()
const shortnames = Dict(
  (".tif", ".tiff") => "GTiff",
  (".nc", ".nc4") => "netCDF",
  (".img",) => "HFA",
  (".xyz",) => "XYZ",
  (".shp",) => "ESRI Shapefile",
  (".geojson",) => "GeoJSON",
  (".fgb",) => "FlatGeobuf",
  (".gdb",) => "OpenFileGDB",
  (".gml",) => "GML",
  (".gpkg",) => "GPKG"
)

## corresponding functions
function find_shortname(fn::AbstractString)
  _, ext = splitext(fn)
  for (k, v) in shortnames
    ext in k && return v
  end
  error("Cannot determine GDAL Driver for $fn")
end

# copied from `GeoArrays`
const gdt_conversion = Dict{DataType,DataType}(
  Bool => UInt8,
  Int8 => UInt8,
  UInt64 => UInt32,
  Int64 => Int32
)

"""Converts type of Array for one that exists in GDAL."""
function cast_to_gdal(A::AbstractArray{<:Real})
  type = eltype(A)
  if type in keys(gdt_conversion)
    newtype = gdt_conversion[type]
    @warn "Casting $type to $newtype to fit in GDAL."
    return newtype, convert(Array{newtype}, A)
  else
    error("Can't cast $(eltype(A)) to GDAL.")
  end
end
