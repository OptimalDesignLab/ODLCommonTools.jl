# a better interface for preallocation and populating sparse matrices

# is SparseMatrixCSC exported from base?
# Ti is type of indices
# Tv is typeof of values
import Base.SparseMatrixCSC
function SparseMatrixCSC{Ti}(sparse_bnds::AbstractArray{Ti, 2}, Tv::DataType)
  # preallocate matrix based on maximum, minimum non zero
  # rows in each column
  # the type of sparse_bnds is used for the indicies
  # the type of val is used for the values
  # the value of val itself is never used

  (tmp, n) = size(sparse_bnds)
  num_nz = 0  # accumulate number of non zero entries

  m = maximum(sparse_bnds)  # get number of rows
  colptr = Array(Ti, n+1)

  if sparse_bnds[1,1] != 0
    colptr[1] = 1
  else
    colptr[1] = 0
  end

  # count number of non zero entries, assign column pointers
  for i=2:(n+1)
    min_row = sparse_bnds[1, i-1]
    max_row = sparse_bnds[2, i-1]

    num_nz += max_row - min_row + 1
    colptr[i] = num_nz + 1
  end

  rowval = zeros(Ti, num_nz)
  nzval = zeros(Tv, num_nz)

  # populate rowvals
  pos = 1
  for i=1:n
    num_vals_i = colptr[i+1] - colptr[i]
    min_row = sparse_bnds[1, i]
#    max_row = sparse_bnds[2, i]

    # write row values to row values
    for j=1:num_vals_i
      rowval[pos] = min_row + j - 1
      pos += 1
    end
  end

  @assert pos == num_nz + 1  # check for sanity

  return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end


import Base.getindex
function getindex{T}(A::SparseMatrixCSC{T}, i::Integer, j::Integer)
# get a nonzero value from A
# for speed, no bounds checking

#  println("using custom getindex")


  row_start = A.colptr[j]
  row_end = A.colptr[j+1] - 1
  row_min = A.rowval[row_start]
  row_max = A.rowval[row_end]

  if i < row_min || i > row_max
    return zero(eltype(A.nzval))
  end

#=
  if i < row_min
    println("in getindex")
   println("i <= row_min, i = ", i, ", j = ", j, "row_min = ", row_min)
   @assert(i >= row_min, string("i = ", i, ", j = ", j, ", row_min = ", row_min))
 end

  if i > row_max
    println("in getindex")
    println("i >= row_max, i = ", i, ", j = ", j, "row_min = ", row_min, ", row_max = ", row_max)
    @assert( i <= row_max)
  end
=#

  offset = i - row_min  # offset due to row

  valindex = row_start + offset

  return A.nzval[valindex]
end

import Base.setindex!
function setindex!{T, Ti}(A::SparseMatrixCSC{T, Ti}, v, i::Integer, j::Integer)
# get a nonzero value from A
# for speed, no bounds checking

#  println("using custom setindex")

  row_start = A.colptr[j]
  row_end = A.colptr[j+1] - 1
  row_min = A.rowval[row_start]
  row_max = A.rowval[row_end]


  if i < row_min || i > row_max
    println(STDERR, "Warning: Cannot change sparsity pattern of this matrix")
    return A
  end
#=  
  if i < row_min
   println("in setindex")
   println("i <= row_min, i = ", i, ", j = ", j, "row_min = ", row_min)
   @assert(i >= row_min)
 end

  if i > row_max
    println("in setindex")
    println("i >= row_max, i = ", i, ", j = ", j, "row_min = ", row_min, ", row_max = ", row_max)
    @assert( i <= row_max)
  end
=#



  offset = i - row_min  # offset due to row

  valindex = row_start + offset
  A.nzval[valindex] = v
  return A
end



  
