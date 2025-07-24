export Point, Line
export line_start, line_end


abstract type AbstractPoint{T} end

Base.@kwdef mutable struct Point3{T} <: AbstractPoint{T}
  x::T
  y::T
  z::T
end

Base.@kwdef mutable struct Point{T} <: AbstractPoint{T}
  x::T
  y::T
end

Base.@kwdef mutable struct Line{T}
  origin::Point{T} = Point(0.0, 0.0)
  azimuth::T = 0.0 # deg, 正北为0, 顺时针为正
  length::T = 2.0 # in deg
  k::T = azimuth2slope(azimuth)
end


@inline line_start(line::Line) = line.origin

function line_end(line::Line)
  p0 = line.origin
  length = line.length
  x, y = p0.x, p0.y
  θ = line.azimuth |> azimuth2deg |> deg2rad
  Point(x + cos(θ) * length, y + sin(θ) * length)
end

function st_bbox(line::Line)
  p0 = line_start(line)
  p1 = line_end(line)
  xmax = max(p0.x, p1.x)
  xmin = min(p0.x, p1.x)
  ymax = max(p0.y, p1.y)
  ymin = min(p0.y, p1.y)
  bbox(; xmin, ymin, xmax, ymax)
end

function is_vertical(line::Line; eps=1e-4)
  abs(mod(line.azimuth, 180)) <= eps # 0或180，认为是垂线
end


rm_empty(xs::AbstractVector) = map(x -> x, filter(!isnothing, xs))


function earth_dist(p1::Point3{T}, p2::Point3{T}) where {T}
  earth_dist((p1.x, p1.y), (p2.x, p2.y))
end
