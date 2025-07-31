export direction, direction_E, direction_N
export azimuth2deg, azimuth2slope, azimuth2xy
export deg2xy


"""
方位角转换为数学角度，限定在0-360°

```bash
  0 ->   90
 90 ->    0
180 ->  -90
270 -> -180
```
"""
function azimuth2deg(Φ::Real)
  Φ = -Φ + 90
  (Φ > 360) && (Φ -= 360)
  (Φ < 0) && (Φ += 360)
  Φ
end

azimuth2slope(ψ) = tan(deg2rad(azimuth2deg(ψ)))

function azimuth2xy(Φ::Real; verbose::Bool=false)
  θ = azimuth2deg(Φ) |> deg2rad
  x, y = cos(θ), sin(θ)
  
  if verbose
    dir = direction(x, y)
    println(dir)
  end
  x, y
end

function deg2xy(θ::Real)
  x, y = cos(θ), sin(θ)
  # dir = direction(x, y)
  # println(dir)
  x, y
end


function direction_E(x::Real; eps=1e-4)
  dir = x > 0 ? "东" : "西"
  abs(x) <= eps ? "正" : dir
end

function direction_N(y::Real; eps=1e-4)
  dir = y > 0 ? "北" : "南"
  abs(y) <= eps ? "正" : dir
end

function direction(x::T, y::T) where {T}
  E = direction_E(x)
  N = direction_N(y)
  E * N
end
