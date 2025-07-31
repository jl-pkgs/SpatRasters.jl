using NearestNeighbors
using LinearAlgebra


"""
计算二次方程系数
"""
function calc_abc(p1::P, p2::P, p3::P, p4::P, out_x::T, out_y::T) where {T <: Real, P <: AbstractPoint{T}}
  # 参考点之间的经向分离
  x_21 = p2.x - p1.x
  x_31 = p3.x - p1.x
  x_42 = p4.x - p2.x

  y_21 = p2.y - p1.y
  y_31 = p3.y - p1.y
  y_42 = p4.y - p2.y

  a = x_31 * y_42 - y_31 * x_42
  b = out_y * (x_42 - x_31) - out_x * (y_42 - y_31) +
      x_31 * p2.y - y_31 * p2.x +
      y_42 * p1.x - x_42 * p1.y
  c = out_y * x_21 - out_x * y_21 + p1.x * p2.y - p2.x * p1.y
  return a, b, c
end


function solve_quadratic(a::T, b::T, c::T; min::T=T(0), max::T=T(1)) where {T<:Real}
  # a x^2 + b x + c = 0
  Δ = b * b - 4 * a * c
  # @show Δ
  # @show max(Δ, T(0))

  sqrt_Δ = sqrt(Base.max(Δ, T(0)))
  # 求解二次多项式
  x_1 = (-b + sqrt_Δ) / (2 * a)
  x_2 = (-b - sqrt_Δ) / (2 * a)
  x_3 = -c / b # (线性情况, a = 0)

  # 找到有效解，即 0 <= t <= 1
  t = _is_out_range(x_1, min, max) ? x_2 : x_1
  t = _is_out_range(x_1, min, max) ? x_3 : t
  t = _is_out_range(x_1, min, max) ? T(NaN) : t
  return t
end


function _is_out_range(x::T, min::T, max::T) where {T}
  isnan(x) && return true
  (x < min || x > max) && return true
end


"""
获取不规则情况的分数距离
"""
function get_fractional_distances_irregular(p1::P, p2::P, p3::P, p4::P, out_x::T, out_y::T) where {T <: Real, P <: AbstractPoint{T}}
  a, b, c = calc_abc(p1, p2, p3, p4, out_y, out_x)
  t = solve_quadratic(a, b, c)
  # 注意p2和p3交换
  s = solve_another_fractional_distance(t, p1.y, p3.y, p2.y, p4.y, out_y)
  return t, s
end


function solve_another_fractional_distance(t::T, y_1::T, y_2::T, y_3::T, y_4::T,
  out_y::T) where T<:Real

  y_21 = y_2 - y_1
  y_43 = y_4 - y_3

  denominator = y_43 * t - y_21 * t + y_3 - y_1
  s = (out_y - y_1 - y_21 * t) / denominator
  @show denominator, s

  # 将值限制在区间[0, 1]，处理除零情况
  if abs(denominator) < eps(T) || _is_out_range(s, T(0), T(1)) 
    s = T(NaN)
  end
  return s
end


export get_fractional_distances_irregular
