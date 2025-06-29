# using D8 method to calculate the flow direction
function D8Direction(dem::AbstractMatrix{T}, i::Int, j::Int, spill::T) where {T<:Real}
  nrow, ncol = size(dem)
  max = 0.0
  steepestSpill = 0
  last_index = 1 # default the first
  grad = 0.0

  for k in 1:8
    i2 = row_goto(i, k)
    j2 = col_goto(j, k)

    # 确保水能够流出去
    !InGrid(i2, j2, nrow, ncol) && continue
    if !isnan(dem[i2, j2]) && (iSpill = dem[i2, j2]) < spill
      grad = (spill - iSpill) / DirLength(k)

      if max < grad
        max = grad
        steepestSpill = k
      end
    end

    if !InGrid(i2, j2, nrow, ncol) || isnan(dem[i2, j2])
      last_index = k # 一般不会发挥作用
    end
  end

  # 这里可能会找不到坡度
  steepestSpill != 0 ? DIR[steepestSpill] : DIR[last_index]
end
