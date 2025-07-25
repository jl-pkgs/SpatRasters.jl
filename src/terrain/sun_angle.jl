include("angle.jl")
using Dates


export angle_SolarDeclination, angle_SunAzimuth, angle_SunElevation
export time2local


utc2local(time::DateTime, lon::Real) = time + dhour(lon / 15)

local2utc(time_local::DateTime, lon::Real) = time_local - dhour(lon / 15)

"北京时间修正到localtime"
function time2local(time::DateTime, lon::Real; lon_ref=120)
  minutes = round(Int, (lon - lon_ref) / 15 * 60) # minutes
  time + Minute(minutes)
end


"""
    angle_SolarDeclination(J; to_deg=false)

# Arguments
- σ: Solar Declination Angle, 黄赤交角（太阳赤纬角）
"""
angle_SolarDeclination(J::Integer) = 0.409 * sin(2pi / 365 * J - 1.39)  # Allen 1998, Eq. 24


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
  E = -cos(σ) * sin(ω) / cos(H) # sinA
  N = (sin(σ) * cos(ψ) - cos(σ) * sin(ψ) * cos(ω)) / cos(H) # cosA

  # 使用 atan2 计算方位角，确保在正确的象限
  A = atan(N, E) # first is y, x
  A = mod(pi/2 - A, 2pi)
  (; H=rad2deg(H), Φ=rad2deg(A))
end


# """
# 任一点，任一时刻的太阳高度角
# """
# function angle_SunElevation(lat::Real, time_local::DateTime)
#   ψ = deg2rad(lat)  # 纬度转弧度
#   J = dayofyear(time_local)
#   σ = angle_SolarDeclination(J)
#   dh = hour(time_local) + minute(time_local) / 60 - 12.0
#   ω = deg2rad(dh * 15) # 时角，上午为负，下午为正

#   sinH = cos(ψ) * cos(σ) * cos(ω) + sin(ψ) * sin(σ)
#   return asin(sinH)
# end
