using Test, SpatRasters, ArchGDAL

import Random: seed!
set_seed(seed) = seed!(seed)

# import NaNStatistics
# using Ipaper: NanQuantile_low, NanQuantile_low!
include("hydro/test-flowdir.jl")

include("test-gdal_polygonize.jl")
include("test-bbox.jl")
include("test-rast.jl")
include("test-st_extract.jl")
include("test-st_mosaic.jl")
include("test-write_gdal.jl")

# include("sf/test_sf.jl")

# using Distributions
# include("Statistics/test-Statistics.jl")

# println(dirname(@__FILE__))
# println(pwd())
# cd(dirname(@__FILE__)) do

## Ipaper
# include("test-agg.jl")
# include("test-par.jl")
# include("test-Ipaper.jl")
# include("test-missing.jl")
# include("test-Pipe.jl")
# include("test-string.jl")
# include("test-list.jl")
# include("test-date.jl")
# include("test-r_base.jl")
# include("test-tools.jl")
