export direction, direction_E, direction_N
export azimuth2deg, azimuth2slope, azimuth2xy
export angle_SolarDeclination, angle_SunAzimuth, angle_SunElevation


using Dates

"""
方位角转换为数学角度，限定在0-360°

```bash
  0 ->   90
 90 ->    0
180 ->  -90
270 -> -180
```
"""
function azimuth2deg(Φ)
  Φ = -Φ + 90
  (Φ > 360) && (Φ -= 360)
  (Φ < 0) && (Φ += 360)
  Φ
end

azimuth2slope(ψ) = tan(deg2rad(azimuth2deg(ψ)))

function azimuth2xy(Φ)
  θ = azimuth2deg(Φ) |> deg2rad
  x, y = cos(θ), sin(θ)

  dir = direction(x, y)
  println(dir)
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


utc2local(time::DateTime, lon::Real) = time + dhour(lon / 15)

local2utc(time_local::DateTime, lon::Real) = time_local - dhour(lon / 15)


"""
    angle_SolarDeclination(J; to_deg=false)

# Arguments
- σ: Solar Declination Angle, 黄赤交角（太阳赤纬角）
"""
angle_SolarDeclination(J::Integer) = 0.409 * sin(2pi / 365 * J - 1.39)  # Allen 1998, Eq. 24


"""
任一点，任一时刻的太阳高度角
"""
function angle_SunElevation(lat::Real, time_local::DateTime)
  ψ = deg2rad(lat)  # 纬度转弧度
  J = dayofyear(time_local)
  σ = angle_SolarDeclination(J)
  dh = hour(time_local) + minute(time_local) / 60 - 12.0
  ω = deg2rad(dh * 15) # 时角，上午为负，下午为正

  sinH = cos(ψ) * cos(σ) * cos(ω) + sin(ψ) * sin(σ)
  return asin(sinH)
end


"""
太阳方位角，正北为0°，正东为90°
"""
function angle_SunAzimuth(lat::Real, time_local::DateTime)
  ψ = deg2rad(lat)  # [deg] -> [radian]
  J = dayofyear(time_local)
  σ = angle_SolarDeclination(J)  # 太阳赤纬角, [radian]

  # 计算时角
  dh = hour(time_local) + minute(time_local) / 60 - 12.0
  ω = deg2rad(dh * 15)  # 时角，上午为负，下午为正

  # 计算太阳高度角
  sinH = cos(ψ) * cos(σ) * cos(ω) + sin(ψ) * sin(σ)
  H = asin(sinH)  # 太阳高度角

  # 计算太阳方位角
  sinA = cos(σ) * sin(ω) / cos(H)
  cosA = (sin(σ) * cos(ψ) - cos(σ) * sin(ψ) * cos(ω)) / cos(H)

  # 使用 atan2 计算方位角，确保在正确的象限
  A = atan(sinA, cosA)
  (; H=rad2deg(H), Φ=rad2deg(A))
end
