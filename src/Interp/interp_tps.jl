using LinearAlgebra
using Base.Threads
export solve_tps, interp_tps, tps_kernel, predict


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
@with_kw struct ThinPlateSpline{TX<:AbstractMatrix,Ty,TΦ}
  λ::Real  # scalar, Stiffness
  x::TX # [n, nd], control points
  y::Ty  # [n, ntime], Homogeneous control point coordinates
  Φ::TΦ  # [n, n], TPS kernel
  d::Ty  # [1+nd, ntime], Affine component
  c::Ty  # [n, ntime], Non-affine component
end

# Thin-plate splines.
# based on http://en.wikipedia.org/wiki/Thin_plate_spline

# tps_basis is applied to the result of a norm which is either positive or zero
# using ifelse is faster than the (? x:y) notation
"""
    Radial Basis Function (径向基函数)

U(r) = r^2 log(r)
"""
function rbf_tps(r::T; eps=1e-6)::T where {T}
  @fastmath r < eps ? T(0) : r * r * log(r)
end

# x: matrix of size KxD
function tps_kernel(x::AbstractMatrix{FT}; distance=distance_norm) where {FT}
  n = size(x, 1) # [n, nx]
  Φ = zeros(FT, n, n)

  @inbounds for i in 1:n              # 上三角
    p1 = @view x[i, :]
    for j in i+1:n
      p2 = @view x[j, :]
      r = distance(p1, p2)
      Φ[i, j] = rbf_tps(r)
    end
  end
  @inbounds for i in 1:n, j in 1:i-1  # 下三角
    Φ[i, j] = Φ[j, i]
  end
  Φ
end

## TODO: 若不是lon, lat，包含了其他协变量，U应该如何计算
function tps_kernel(x1::AbstractMatrix{FT}, x2::AbstractMatrix{FT}; distance=distance_norm) where {FT}
  n_control = size(x1, 1)
  n_points = size(x2, 1)
  Φ = zeros(FT, n_points, n_control)

  @inbounds for i in 1:n_points
    p2 = @view x2[i, :]
    for j in 1:n_control
      p1 = @view x1[j, :]
      r = distance(p1, p2)
      Φ[i, j] = rbf_tps(r)
    end
  end
  Φ
end

"""
	solve_tps(x,y,λ,compute_affine=true)

find solution of tps transformation 
(required for some operations but takes additional time.)

## Arguments
	`x`: control points, n
	`y`: deformed (warped) control points
	`λ`: stiffness coefficient

returns a `ThinPlateSpline` structure which can be supplied as an argument to `tps_deform`.

## References
1. Ghosh A, Kindermann D D S. Efficient thin plate spline interpolation and its
   application to adaptive optics[M]. na, 2010.

## See almost
	`ThinPlateSpline`
"""
function solve_tps(x::AbstractMatrix, y::AbstractArray, λ::Real; distance::Function=distance_norm)
  n, D = size(x)

  # homogeneous coordinates
  X = hcat(ones(n, 1), x)
  Y = y
  Φ = tps_kernel(x; distance) # compute TPS kernel, also known as K

  # full QR decomposition
  Q, r = qr(X)
  Q1 = Q[:, 1:(D+1)]
  Q2 = Q[:, (D+2):end]

  # warping coefficients
  b = inv(UniformScaling(λ) + Q2' * Φ * Q2) * Q2' * Y # Ghosh 2010, Eq. 3.7
  c = Q2 * b
  d = r \ (Q1' * (Y - Φ * c))                         # Ghosh 2010, Eq. 3.6
  return ThinPlateSpline(λ, x, y, Φ, d, c)
end


function predict(tps::ThinPlateSpline, x2::AbstractMatrix{FT};
  distance::Function=distance_norm, progress=true) where {FT}
  (; d, c) = tps
  x1 = tps.x

  n_points = size(x2, 1)   # 输入点数
  n_control = size(x1, 1)  # 控制点数
  nx = size(x1, 2)         # 空间维度
  ntime = size(tps.y, 2)   # 输出维度

  R = zeros(FT, n_points, ntime) # 输出初始化
  U = tps_kernel(x1, x2; distance) # [n_points, n_control]

  p = Progress(ntime)
  @inbounds for j in 1:ntime
    progress && next!(p)

     @threads for i in 1:n_points
      z = d[1, j] # 常数项
      for l in 1:nx
        z += d[l+1, j] * x2[i, l]
      end

      for ctrl in 1:n_control
        # _x1 = @view x1[ctrl, :]
        # _x2 = @view x2[i, :]
        # r = distance(_x1, _x2) # 可以提前算好mat_U
        z += U[i, ctrl] * c[ctrl, j]
      end
      R[i, j] = z
    end
  end
  return R
end


function interp_tps(X::AbstractMatrix, Y::AbstractArray, target::SpatRaster;
  λ=0.01, distance::Function=distance_norm, progress=true, kw...)

  ntime = size(Y, 2)
  tps = solve_tps(X, Y, λ)
  
  lon, lat = st_dims(target)
  Lon, Lat = meshgrid(lon, lat)
  
  X2 = [Lon[:] Lat[:]]
  nlon, nlat, _ = size(target)
  R = predict(tps, X2; distance, progress, kw...)
  rast(reshape(R, nlon, nlat, ntime), target)
end
