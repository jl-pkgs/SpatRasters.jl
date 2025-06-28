abstract type AbstractSpatRaster{T,N} end

Base.@kwdef mutable struct SpatRaster{T,N} <: AbstractSpatRaster{T,N}
  A::AbstractArray{T,N}
  b::bbox = bbox(-180.0, -90.0, 180.0, 90.0)
  cellsize::NTuple{2,Real}
  lon::AbstractVector{<:Real}
  lat::AbstractVector{<:Real}
  time::Union{AbstractVector,Nothing} = nothing
  bands::Union{AbstractVector{String},Nothing} = nothing
  name::String = "Raster"
  nodata::Union{AbstractVector{T},Nothing} = nothing
end

function nband(ra::SpatRaster{T,N}) where {T,N}
  n = 1
  N == 3 && (n == size(ra.A)[end])
  N > 3 && (n = size(ra.A)[end])
  return n
end

"""
    SpatRaster(A, b::bbox; reverse_lat=true, kw...)

- `kw`: other parameters: `time`, `name`
"""
function SpatRaster(A::AbstractArray{T,N}, b::bbox=bbox(-180.0, -90.0, 180.0, 90.0);
  reverse_lat=true, time=nothing, bands=nothing, name="Raster", kw...) where {T,N}

  if N == 3 && size(A, 3) == 1
    A = A[:, :, 1]
  end
  cellsize = bbox2cellsize(b, size(A))
  lon, lat = bbox2dims(b; cellsize, reverse_lat)
  SpatRaster(; A, b, cellsize, lon, lat, time, bands, name, kw...)
end

function SpatRaster(A::AbstractArray, r::SpatRaster; reverse_lat=true)
  (; b, time, bands, name) = r
  if size(A)[1:2] != size(r.A)[1:2]
    lon, lat = bbox2dims(b; size=size(A), reverse_lat)
    cellsize = bbox2cellsize(b, size(A))
  else
    (; lon, lat, cellsize) = r
  end
  SpatRaster(; A, b, cellsize, lon, lat, time, bands, name) # rebuild
end

function SpatRaster(f::String; kw...)
  A = read_gdal(f)
  bands = bandnames(f)
  nodata = gdal_nodata(f)
  SpatRaster(A, st_bbox(f); bands, nodata, kw...)
end


function Base.getindex(ra::AbstractSpatRaster, i, j, args...; deep=true)
  (; A, cellsize, lon, lat, time, bands, name, nodata) = ra

  cols = repeat([:], max(ndims(ra) - 2 - length(args), 0))
  inds = (i, j, args..., cols...)
  
  if length(args) > 0
    k = args[1]
    !isnothing(bands) && (bands = bands[k])
    !isnothing(time) && (time = time[k])
  end

  lon, lat = st_dims(ra)
  _lon, _lat = lon[i], lat[j]
  _b = st_bbox(_lon, _lat)
  _A = deep ? A[inds...] : @view A[inds...]

  rast(; A=_A, b=_b, lon=_lon, lat=_lat, time, bands,
    cellsize, name, nodata)
end

function Base.getindex(ra::AbstractSpatRaster, i::Int, j::Int, args...)
  cols = repeat([:], max(ndims(ra) - 2 - length(args), 0))
  inds = (i, j, args..., cols...)
  ra.A[inds...]
end


Base.ndims(ra::AbstractSpatRaster) = ndims(ra.A)
Base.size(ra::AbstractSpatRaster) = size(ra.A)
Base.size(ra::AbstractSpatRaster{T,2}) where {T} = (size(ra.A)..., 1)
# Base.parent(ra::AbstractSpatRaster) = ra.A
# Base.iterate(ra::AbstractSpatRaster) = iterate(ra.A)
# Base.length(ra::AbstractSpatRaster) = length(ra.A)
# Base.size(ra::AbstractSpatRaster) = size(ra.A)
# Base.eltype(::Type{AbstractSpatRaster{T}}) where {T} = T
# Base.map(f, ra::AbstractSpatRaster) = SpatRaster(map(f, ra.A), ra)

# !note about NaN values
Base_ops = ((:Base, :+), (:Base, :-), (:Base, :*), (:Base, :/),
  (:Base, :>), (:Base, :<), (:Base, :>=), (:Base, :<=),
  (:Base, :!=),
  (:Base, :&), (:Base, :|))

for (m, f) in Base_ops
  # _f = Symbol(m, ".:", f)
  @eval begin
    $m.$f(a::AbstractSpatRaster, b::AbstractSpatRaster) = begin
      size(a) != size(b) && throw(DimensionMismatch("size mismatch"))
      SpatRaster($m.$f.(a.A, b.A), a)
    end

    $m.$f(a::AbstractSpatRaster, b::Real) = SpatRaster($m.$f.(a.A, b), a)
    $m.$f(a::Real, b::AbstractSpatRaster) = SpatRaster($m.$f.(a, b.A), b)
  end
end


import Base: ==
function ==(x::SpatRaster, y::SpatRaster)
  x.b == y.b && x.A == y.A && x.nodata == y.nodata
end


function Base.show(io::IO, x::SpatRaster)
  T = eltype(x.A)
  printstyled(io, "SpatRaster{$T}: ", color=:blue)
  printstyled(io, "$(x.name)\n", color=:green, underline=true)

  print(io, "  A        : ")
  obj_size(x.A)

  println(io, "  b        : $(x.b)")
  println(io, "  cellsize : $(x.cellsize)")
  println(io, "  lon, lat : $(x.lon), $(x.lat)")
  if !isnothing(x.time)
    time_beg = x.time[1]
    time_end = x.time[end]
    println(io, "  time     : $time_beg ~ $time_end, ntime=$(length(x.time))")
  end
  println(io, "  bands    : $(x.bands)")
  print(io, "  nodata   : $(x.nodata)")
  return nothing
end

rast = SpatRaster

export AbstractSpatRaster, SpatRaster, rast
