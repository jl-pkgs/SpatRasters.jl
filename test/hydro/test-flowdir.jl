using SpatRasters, ArchGDAL, Test

@testset "flowdir" begin
  dem = read_gdal(guanshan_dem, 1)
  @time dir_julia = FillDEM_FlowDirection(dem) |> gis2tau

  dir_cpp = read_gdal(guanshan_flowdir_cpp, 1) |> gis2tau
  @test dir_cpp == dir_julia
end

# b = st_bbox(f)
# lon, lat = st_dims(f)
# dem_bak = deepcopy(dem)
# obj_size(dem)
