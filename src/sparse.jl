# Description:  
#   A better interface for preallocation and populating sparse matrices

# is SparseMatrixCSC exported from base?
# Ti is type of indices
# Tv is typeof of values
import Base.SparseMatrixCSC

function SparseMatrixCSC{Tv}(mesh::AbstractMesh, ::Type{Tv})
# construct a SparseMatrixCSC that preallocates the space needed for values
# tightly
# this should be exact for DG and a slight overestimate for CG


  # calculate lengths of vectors
  # nel_per_el: number of elements (including self) each element is connected to
  if mesh.coloringDistance == 0
    nel_per_el = 1
  elseif mesh.coloringDistance == 2
    nel_per_el = size(mesh.neighbor_nums, 1)
  else
    throw(ErrorException("Unsupported coloring distance"))
  end
  ndof = mesh.numDof
  nDofPerElement = mesh.numNodesPerElement*mesh.numDofPerNode
  nvals_per_column = mesh.numNodesPerElement*mesh.numDofPerNode*nel_per_el
  nvals = nvals_per_column*mesh.numDof

  colptr = Array(Int64, ndof+1)
  rowvals = zeros(Int64, nvals)
  nzvals = zeros(Tv, nvals)

  # set up colptr
  colptr[1] = 1
  for i=2:length(colptr)
    colptr[i] = colptr[i-1] + nvals_per_column
  end

  # set up rowvals
  dofs_i = zeros(eltype(mesh.dofs), nvals_per_column)
  elnums_i = zeros(eltype(mesh.pertNeighborEls), nel_per_el)
  # loop over elements because all nodes on an element have same sparsity
  for i=1:mesh.numEl

    # get the element numbers
    pos = 1
    for j=1:size(mesh.pertNeighborEls, 2)
      val = mesh.pertNeighborEls[i, j]
      if val > 0
        elnums_i[pos] = val
        pos += 1
      end
    end

    # copy dof numbers into array
    for j = 1:(pos-1)
      el_j = elnums_i[j]
      dofs_j = view(mesh.dofs, :, :, el_j)
      src_j = view(dofs_i, ((j-1)*nDofPerElement+1):(j*nDofPerElement))
      copyDofs(dofs_j, src_j)
    end

    dofs_used = view(dofs_i, 1:(pos-1)*nDofPerElement)
    # sort them
    sort!(dofs_used)
    @assert dofs_used[1] != 0

    # note: although it would be technically easy to de-duplicate the dof
    #       numbers here, we don't because that would make it difficult to 
    #       calculate the starting position on rowvals for each dof

    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        # compute starting location in rowvals
        mydof = mesh.dofs[k, j, i]
        start_idx = (mydof-1)*nvals_per_column + 1
        # set them in rowvals
        for p=1:length(dofs_used)
          rowvals[start_idx + p - 1] = dofs_used[p]
        end
      end
    end

  end  # end loop over elements

  return SparseMatrixCSC(ndof, ndof, colptr, rowvals, nzvals)

end

function copyDofs{T}(src::AbstractArray{T, 2}, dest::AbstractArray{T, 1})
  pos = 1
  for i=1:size(src, 2)
    for j=1:size(src, 1)
      dest[pos] = src[j, i]
      pos += 1
    end
  end
end

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

const band_dense = false

if band_dense
  # setindex for dense within the band matrix
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


else
  function setindex!{T, Ti}(A::SparseMatrixCSC{T, Ti}, v, i::Integer, j::Integer)
    row_start = A.colptr[j]
    row_end = A.colptr[j+1] - 1
    rowvals_extract = unsafe_view(A.rowval, row_start:row_end)
    val_idx = fastfind(rowvals_extract, i)
    idx = row_start + val_idx - 1
    A.nzval[idx] = v

    return A


  end

  function getindex{T}(A::SparseMatrixCSC{T}, i::Integer, j::Integer)
    row_start = A.colptr[j]
    row_end = A.colptr[j+1] - 1
    rowvals_extract = unsafe_view(A.rowval, row_start:row_end)
    val_idx = fastfind(rowvals_extract, i)
    idx = row_start + val_idx -1
    return A.nzval[idx]
   
  end

end  # end if band_dense
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


