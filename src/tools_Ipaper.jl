export isscalar

isscalar(x) = !isa(x, AbstractArray)
isscalar(::Nothing) = false


file_ext(file::String) = file[findlast(==('.'), file):end]


"""
    obj_size(x)
    obj_size(dims, T)

# Examples
```julia
dims = (100, 100, 200)
x = zeros(Float32, dims)
obj_size(x)
obj_size(dims, Float32)
```
"""
function obj_size(io::IO, x)
  ans = Base.summarysize(x) / 1024^2
  ans = round(ans, digits=2)
  print(io, typeof(x), " | ", size(x), " | ")
  printstyled(io, "$ans Mb\n"; color=:blue, bold=true, underline=true)
end

function obj_size(dims, T)
  ans = Base.summarysize(T(0)) * prod(dims) / 1024^2
  ans = round(ans, digits=2)
  print(T, " | ", dims, " | ")
  printstyled("$ans Mb\n"; color=:blue, bold=true, underline=true)
end


# findnear(x::Real, vals::AbstractVector) = argmin(abs.(vals .- x))
# function findnear(x::Real, y::Real, lon::AbstractVector, lat::AbstractVector)
#   i = findnear(x, lon)
#   j = findnear(y, lat)
#   return i, j
# end
function findnear(x::Real, vals::AbstractVector; cell::Real=NaN, tol=1e-2)
  diff = abs.(vals .- x)
  i = argmin(diff)
  isnan(cell) && return i
  diff[i] <= (0.5 + tol) * abs(cell) ? i : -1 # 在1个网格内
end

# cellsize需要是准确的
function findnear((x, y)::Tuple{Real,Real}, lon::AbstractVector, lat::AbstractVector;
  cellx::Real=NaN, celly::Real=NaN, tol=1e-2)
  i = findnear(x, lon; cell=cellx, tol)
  j = findnear(y, lat; cell=celly, tol)
  (i == -1 || j == -1) && (return nothing)
  return i, j
end

findnear(x::Real, y::Real, lon::AbstractVector, lat::AbstractVector; cellx::Real=NaN, celly::Real=NaN, tol=1e-2) =
  findnear((x, y), lon, lat; cellx, celly, tol)
