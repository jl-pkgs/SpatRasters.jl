using Ipaper
using Ipaper.sf
using Test
using ArchGDAL
import GeoDataFrames as GDF

@testset "gdal_polygonize" begin
  set_seed(1)
  n = 2
  A = rand(Bool, 7 * n, 4 * n)
  b = bbox([70, 15, 140, 55]...)

  A[1] = 0
  ra = rast(A, b)
  ra.A

  f = "test_nodata.tif"
  # st_write(ra, f; nodata=0);
  write_gdal(ra, f) # both works

  # gdal_info(f)
  gdal_polygonize(f, "out.shp", mask=true, diag=false)
  df = GDF.read("out.shp")
  @test size(df, 1) == 15

  gdal_polygonize(f, "out.shp", mask=true, diag=true)
  df = GDF.read("out.shp")
  @test size(df, 1) == 3

  gdal_polygonize(f, "out.shp", mask=false, diag=false)
  df = GDF.read("out.shp")
  @test size(df, 1) == 31

  gdal_polygonize(f, "out.shp", mask=false, diag=true)
  df = GDF.read("out.shp")
  @test size(df, 1) == 9
  
  rm_shp("out.shp")
  rm(f)
end


@testset "gdal_polygonize multibands" begin
  set_seed(1)
  n = 10
  ntime = 100
  A = rand(Bool, 7 * n, 4 * n, 100)
  b = bbox([70, 15, 140, 55]...)

  ra = rast(A, b)
  ra.A

  ## 1. without nodata
  f = "test.tif"
  st_write(ra, f) # both works

  fout = "out.gdb"
  @time gdal_polygonize(f, fout, mask=true, diag=false)
  @test nlayer(fout) == 100
  ogr_info(fout)
  run(`rm -rf $fout`)

  fout = "out.gdb"
  @time gdal_polygonize(f, fout; bands=1:2, mask=true, diag=false)
  @test nlayer(fout) == 2
  run(`rm -rf $fout`)

  fout = "out.gpkg" # about 10 times slower
  @time gdal_polygonize(f, fout; bands=1:2, mask=true, diag=false)
  @test nlayer(fout) == 2
  ogr_info(fout)
  run(`rm -rf $fout`)
  rm(f)
end
