function get_four_closest_corners(
  in_x::Vector{T}, in_y::Vector{T}, out_x::Vector{T}, out_y::Vector{T},
  neighbours::Int, index_array::Matrix{Int32}
) where T<:Real

  stride, valid_corners = get_stride_and_valid_corner_indices(
    out_x, out_y, in_x, in_y, neighbours
  )

  res = Vector{Matrix{T}}()
  indices = Vector{Vector{Int32}}()

  for valid in valid_corners
    x__, y__, idx = get_corner(stride, valid, in_x, in_y, index_array)
    push!(res, hcat(x__, y__))
    push!(indices, idx)
  end
  return res, hcat(indices...)
end


"获取步长和有效角点索引"
function get_stride_and_valid_corner_indices(
  out_x::Vector{T}, out_y::Vector{T},
  in_x::Vector{T}, in_y::Vector{T},
  neighbours::Int
) where T<:Real

  out_x_tile = repeat(out_x, 1, neighbours)
  out_y_tile = repeat(out_y, 1, neighbours)

  # 获取两个方向的差异
  x_diff = out_x_tile .- in_x
  y_diff = out_y_tile .- in_y

  stride = 1:size(x_diff, 1) ## 划分4个角
  valid_corners = (
    (x_diff .> 0) .& (y_diff .< 0),  # 左上角
    (x_diff .< 0) .& (y_diff .< 0),  # 右上角  
    (x_diff .> 0) .& (y_diff .> 0),  # 左下角
    (x_diff .< 0) .& (y_diff .> 0)   # 右下角
  )
  return stride, valid_corners
end
