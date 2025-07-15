using LinearAlgebra
using Base.Threads
export ThinPlateSpline, solve_tps, predict


"""
	ThinPlateSpline{Λ,TX<:AbstractMatrix,M,Q}

the type (structure) holding the deformation information. This is needed to apply the deformation to other points using `tps_deform`.
(see http://en.wikipedia.org/wiki/Thin_plate_spline for more information)

# Members
`λ::Λ`  # Stiffness.
`x1::X` # control points
`Y::M`  # Homogeneous control point coordinates
`Φ::M`  # TPS kernel
`d::M`  # Affine component
`c::M`  # Non-affine component
"""
@with_kw struct ThinPlateSpline{Λ,TX<:AbstractMatrix,M,Q}
  λ::Λ  # scalar, Stiffness
  x::TX # [n, nd], control points
  y::M  # [n, ntime], Homogeneous control point coordinates
  Φ::Q  # [n, n], TPS kernel
  d::M  # [1+nd, ntime], Affine component
  c::M  # [n, ntime], Non-affine component
end

# Thin-plate splines.
# based on http://en.wikipedia.org/wiki/Thin_plate_spline

# tps_basis is applied to the result of a norm which is either positive or zero
# using ifelse is faster than the (? x:y) notation
"""
U(r) = r^2 log(r)
"""
U(r::T) where {T} = ifelse(r < eps(r), zero(T), r * r * log(r))

# _norm(a) = sqrt(sum(a .^ 2))

# x: matrix of size KxD
function tps_kernel(x::AbstractMatrix{FT}; distance=distance_norm) where {FT}
  n = size(x, 1) # [n, nx]
  Φ = zeros(FT, n, n)

  @inbounds for i in 1:n              # 上三角
    p1 = @view x[i, :]
    for j in i+1:n
      p2 = @view x[j, :]
      r = distance(p1, p2)
      Φ[i, j] = U(r)
    end
  end
  @inbounds for i in 1:n, j in 1:i-1  # 下三角
    Φ[i, j] = Φ[j, i]
  end
  Φ
end

"""
	tps_solve(x,y,λ,compute_affine=true)

find solution of tps transformation 
(required for some operations but takes additional time.)

# Arguments
	`x`: control points, n
	`y`: deformed (warped) control points
	`λ`: stiffness coefficient
	`compute_affine`: computes affine component if `true`

returns a `ThinPlateSpline` structure which can be supplied as an argument to `tps_deform`.

# See almost
	`ThinPlateSpline`
"""
function solve_tps(x::AbstractMatrix, y::AbstractArray, λ::Real; distance::Function=distance_norm)
  n, D = size(x)

  # homogeneous coordinates
  X = hcat(ones(n, 1), x)
  Y = y
  Φ = tps_kernel(x; distance) # compute TPS kernel

  # full QR decomposition
  Q, r = qr(X)
  q1 = Q[:, 1:(D+1)]
  q2 = Q[:, (D+2):end]

  # warping coefficients
  b = inv(UniformScaling(λ) + q2' * Φ * q2) * q2' * Y # Ghosh 2010, Eq. 3.6
  c = q2 * b
  d = r \ (q1' * (Y - Φ * c))                         # Ghosh 2010, Eq. 3.7
  return ThinPlateSpline(λ, x, y, Φ, d, c)
end


function predict(tps::ThinPlateSpline, x2::AbstractMatrix{FT};
  distance::Function=distance_norm) where {FT}
  (; d, c) = tps
  x1 = tps.x

  n_points = size(x2, 1)   # 输入点数
  n_control = size(x1, 1)  # 控制点数
  nx = size(x1, 2)         # 空间维度
  ntime = size(tps.y, 2)   # 输出维度

  R = zeros(FT, n_points, ntime) # 输出初始化

  @inbounds for j in 1:ntime
    @threads for i in 1:n_points
      z = d[1, j] # 常数项
      for l in 1:nx
        z += d[l+1, j] * x2[i, l]
      end

      for ctrl in 1:n_control
        _x1 = @view x1[ctrl, :]
        _x2 = @view x2[i, :]
        r = distance(_x1, _x2)
        z += U(r) * c[ctrl, j]
      end
      R[i, j] = z
    end
  end
  return R
end
