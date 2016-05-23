# Description:  
#   A better interface for preallocation and populating sparse matrices

# is SparseMatrixCSC exported from base?
# Ti is type of indices
# Tv is typeof of values
import Base.SparseMatrixCSC

function SparseMatrixCSC{Ti}(sparse_bnds::AbstractArray{Ti, 2}, Tv::DataType)
# TODO: @doc this
# preallocate matrix based on maximum, minimum non zero
# rows in each column
# the type of sparse_bnds is used for the indicies
# the type of val is used for the values
# the value of val itself is never used

  println("creating SparseMatrixCSC")

  (tmp, n) = size(sparse_bnds)
  num_nz = 0  # accumulate number of non zero entries

  m = maximum(sparse_bnds)  # get number of rows
  colptr = Array(Int64, n+1)  # should be Ti

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

  rowval = zeros(Int64, num_nz)  # should be Ti
  nzval = zeros(Tv, num_nz)

  # populate rowvals
  pos = 1
  for i=1:n
    num_vals_i = colptr[i+1] - colptr[i]
    min_row = sparse_bnds[1, i]

    # write row values to row values
    for j=1:num_vals_i
      rowval[pos] = min_row + j - 1
      pos += 1
    end
  end

  @assert pos == num_nz + 1  # check for sanity

  println("average bandwidth = ", pos/m)
  return SparseMatrixCSC(m, n, colptr, rowval, nzval)
end

#------------------------------------------------------------------------------
# Access methods
import Base.getindex
import Base.setindex!
import Base.fill!

function getindex{T}(A::SparseMatrixCSC{T}, i::Integer, j::Integer)
# TODO: @doc this
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

  offset = i - row_min  # offset due to row
  valindex = row_start + offset

  return A.nzval[valindex]

end

function setindex!{T, Ti}(A::SparseMatrixCSC{T, Ti}, v, i::Integer, j::Integer)
# TODO: @doc this
# get a nonzero value from A
# for speed, no bounds checking

#  println("using custom setindex")

  row_start = A.colptr[j]
  row_end = A.colptr[j+1] - 1
  row_min = A.rowval[row_start]
  row_max = A.rowval[row_end]

  if i < row_min || i > row_max
    println(STDERR, "Warning: Cannot change sparsity pattern of this matrix")
    println(STDERR, "    i = ", i, ", j = ", j, " value = ", v)
    return A
  end

  offset = i - row_min  # offset due to row
  valindex = row_start + offset
  A.nzval[valindex] = v

  return A

end

function fill!(A::SparseMatrixCSC, val)
  fill!(A.nzval, val)
  return nothing
end

@doc """
### ODLCommonTools.findfast

  This function searches a sorted array for a given value, returning 0
  if it is not found.  

  The algorithm is nearly branchless and performs well compared to
  standard implementations.  

  The search takes a maximum of log2(n) + 2 iterations when the requested
  value is present and n iteration if it is not found.

  Inputs:
    arr: array of integers
    val: value to find

  Outputs:
    idx: the index of the array containing the value, 0 if not found
"""->
function fastfind{T <: Integer}(a::AbstractArray{T}, val)

  foundflag = false
  lbound = 1
  ubound = length(a)
  idx = lbound + div(ubound - lbound, 2)
#  itermax = floor(log2(length(a))) + 2
  itermax = length(a)
  itr = 0


#  println("lbound = ", lbound)
#  println("ubound = ", ubound)

  while ( a[idx] != val && itr <= itermax)
#    println("\ntop of loop, idx = ", idx)
    if a[idx] > val  # the value lies in the left half 
      ubound = idx
      idx = lbound + fld(ubound - lbound, 2)
#      println("updating ubound = ", ubound)
    else  # a[idx] < val  # value lies in the right half
      lbound = idx
      idx = lbound + cld(ubound - lbound, 2)
#      println("updating lbound = ", lbound)
    end
    
    itr += 1
  end

    successflag = (itr <= itermax)
  return idx*successflag
#  return idx
end


