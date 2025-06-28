# using Test, Ipaper, Ipaper.sf, ArchGDAL

# not perfect
function black2white(A::Array{UInt8,3}; threshold::Int=10, white::Int=255)
  threshold = UInt8(threshold)
  white = UInt8(white)
  A = copy(A)
  n, m, p = size(A)
  @inbounds for i = 1:n, j = 1:m
    r = A[i, j, 1]
    g = A[i, j, 2]
    b = A[i, j, 3]

    if r <= threshold && g <= threshold && b <= threshold
      A[i, j, 1] = white
      A[i, j, 2] = white
      A[i, j, 3] = white
    end
  end
  A
end

@testset "write_gdal RGB_tiff" begin
  b = bbox(-180.0, -60.0, 180.0, 90.0)
  A = rand(UInt8, 4, 4, 3)
  r = rast(A, b)

  ## write_tiff, support RGB directly now
  f_data = "test-data.tif"
  f_rast = "test-rast.tif"

  write_gdal(A, f_data)
  write_gdal(r, f_rast)

  gdal_info(f_data)
  @test gdalinfo(f_rast)["bbox"] == b

  isfile(f_data) && rm(f_data)
  isfile(f_rast) && rm(f_rast)
end

# # using Plots
# r = C[:, :, 1][:]
# g = C[:, :, 2][:]
# b = C[:, :, 3][:]

# xlim = (0, 50)
# plot(
#   histogram(r, nbins=100, title="R", xlims=xlim), 
#   histogram(g, nbins=100, title="G", xlims=xlim),
#   histogram(b, nbins=100, title="B", xlims=xlim)
# )
# ## 需要训练一个二分类的算法

# # C2 = black2white(C; threshold=40)
# # write_gdal(C, "test-raw.tif")
# # write_gdal(C2, "test-image.tif")
