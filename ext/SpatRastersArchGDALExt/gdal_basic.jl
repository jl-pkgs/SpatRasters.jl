# nlayer is for org
function nlayer(f)
  ArchGDAL.read(f) do ds
    ArchGDAL.nlayer(ds)
  end
end

function nband(f)
  gdal_open(f) do ds
    GDAL.gdalgetrastercount(ds)
  end
end

# nraster = nband
# nlyr = nband

gdal_close(ds::Ptr{Nothing}) = GDAL.gdalclose(ds)

# gdal_open(file::AbstractString) = ArchGDAL.read(file)
function gdal_open(f::AbstractString, mode=GDAL.GA_ReadOnly, args...)
  GDAL.gdalopen(f, mode, args...)
end

function gdal_open(func::Function, args...; kwargs...)
  ds = gdal_open(args...; kwargs...)
  try
    func(ds)
  finally
    gdal_close(ds)
  end
end

function bandnames(f)
  n = nband(f)
  gdal_open(f) do ds
    map(iband -> begin
        band = GDAL.gdalgetrasterband(ds, iband)
        GDAL.gdalgetdescription(band)
      end, 1:n)
  end
end

# works
function set_bandnames(f, bandnames)
  n = nband(f)
  gdal_open(f, GDAL.GA_Update) do ds
    map(iband -> begin
        band = GDAL.gdalgetrasterband(ds, iband)
        GDAL.gdalsetdescription(band, bandnames[iband])
      end, 1:n)
  end
  return nothing
end

# gdal_nodata(band::ArchGDAL.GDALRasterBand) = ArchGDAL.GDAL.gdalgetrasternodatavalue(band.ptr)
function gdal_nodata(f)
  n = nband(f)
  gdal_open(f, GDAL.GA_Update) do ds
    map(iband -> begin
        band = GDAL.gdalgetrasterband(ds, iband)
        nodata = GDAL.gdalgetrasternodatavalue(band, Ref(Cint(0)))
        Type = ArchGDAL.pixeltype(band)
        Type(nodata)
      end, 1:n)
  end
end

# function gdal_nodata(f::String, band::Int)
#   ArchGDAL.read(f) do dataset
#     band = ArchGDAL.getband(dataset, band)
#     ArchGDAL.getnodatavalue(band)
#   end
# end

function gdal_info(f)
  run(`$(gdalinfo_path()) $f`)
  return nothing
end

function ogr_info(f)
  run(`$(ogrinfo_path()) $f`)
  return nothing
end
