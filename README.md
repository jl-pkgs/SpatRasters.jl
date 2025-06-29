# SpatRasters.jl: Simple Spatial Rasters in Julia

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jl-pkgs.github.io/SpatRasters.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/SpatRasters.jl/dev)
[![CI](https://github.com/jl-pkgs/SpatRasters.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/SpatRasters.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/jl-pkgs/SpatRasters.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jl-pkgs/SpatRasters.jl)

> Dongdong Kong

`SpatRaster` is a simple spatial raster with WGS84 projection, abbreviated as
`rast`.

Due to its simplicity, `SpatRasters.jl` doesn't rely on any other geospatial
packages, which makes `SpatRasters.jl` more lightweight than `GeoArrays` and
`Rasters.jl`.

`SpatRasters.jl` used the class and function of R language `terra` and `sf`
packages.
