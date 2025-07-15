using Test, SpatRasters, ArchGDAL

import Random: seed!
set_seed(seed) = seed!(seed)

# include("hydro/test-flowdir.jl")
include("Interp/test-tps.jl")
include("Interp/test-angle.jl")

include("test-gdal_polygonize.jl")
include("test-bbox.jl")
include("test-rast.jl")
include("test-st_extract.jl")
include("test-st_mosaic.jl")
include("test-write_gdal.jl")

# println(dirname(@__FILE__))
# println(pwd())
# cd(dirname(@__FILE__)) do
