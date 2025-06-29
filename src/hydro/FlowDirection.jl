"""
  fillDEMAndComputeFlowDirection(dem::AbstractMatrix)

#! dem was modified!

fill_depressions_wang_and_liu

Wang&Liu(2006) to fill depressions
"""
function FillDEM_FlowDirection(dem::AbstractMatrix{T}; nchunk = 50) where {T}
  dem = deepcopy(dem)
  nrow, ncol = size(dem)
  isProcessed = falses(nrow, ncol)
  flowdir = zeros(UInt8, nrow, ncol)

  queue = PriorityQueue{Tuple{Int,Int},T}() # Base.Order.Reverse
  n_valid = 0

  # Push all border cells into the queue
  for i in 1:nrow, j in 1:ncol
    isnan(dem[i, j]) && continue
    n_valid += 1

    for k in 1:8
      i2, j2 = row_goto(i, k), col_goto(j, k)
      if !InGrid(i2, j2, nrow, ncol) || isnan(dem[i2, j2])
        queue[i, j] = dem[i, j]
        isProcessed[i, j] = true
        break # ? not bug, 为何跳出for loop
      end
    end
  end

  count = 0
  _spill = 0.0
  
  chunk = n_valid ÷ nchunk
  p = Progress(Int(nchunk))
  
  while !isempty(queue)
    mod(count, chunk) == 0 && next!(p)
    count += 1
    (i, j), spill = peek(queue)
    dequeue!(queue)

    for k in 1:8
      i2, j2 = row_goto(i, k), col_goto(j, k)
      
      !InGrid(i2, j2, nrow, ncol) && continue
      if !isnan(dem[i2, j2]) && !isProcessed[i2, j2]

        iSpill = dem[i2, j2] # next 低地
        if iSpill <= spill
          _spill = spill
          flowdir[i2, j2] = DIR_INV[k] # 站在洼地的角度
        else
          _spill = iSpill # 碰到了洼地
        end

        isProcessed[i2, j2] = true
        dem[i2, j2] = _spill
        queue[i2, j2] = _spill
      end
    end
    
    if (flowdir[i, j] == NODATA)
      flowdir[i, j] = D8Direction(dem, i, j, spill)
    end
  end
  flowdir
end


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


export FillDEM_FlowDirection, D8Direction
