# Description: miscellaneous functions

using ArrayViews

export rmfile, printbacktrace, smallmatvec!, smallmatvec, smallmatTvec!, 
        smallmatTvec, smallmatmat!, 
        smallmatmat, smallmatmatT!, smallmatmatT, checkZeroRows, 
        checkZeroColumns, checkIdenticalColumns, checkSparseColumns,
        checkSparseRows, findLarge, isSymmetric, make_symmetric!

@doc """
 ### Tools rmfile

 Removes the file if it exists.  This is only necessary because it is an error
   to remove a file that does not exist.  Does not work on directories

   Arguments:
     fname::AbstractString  :  name of file

"""->
function rmfile(fname::AbstractString; recursive=false)
  if isfile(fname)
    rm(fname, recursive=recursive)
  end
  return nothing
end

# this function is both slow and dangerous
function readdlmc(fname)
  map(x->eval(parse(x)),readcsv(fname))
end

function printbacktrace()
# print the current call stack
  bt = backtrace()
  s = sprint(io->Base.show_backtrace(io, bt))
  println("backtrace: ", s)

  return s
end


function smallmatvec!{T, T2, T3}(A::AbstractArray{T,2}, 
           x::AbstractArray{T2,1}, b::AbstractArray{T3, 1})
# TODO: this comment needs to be a @doc
# performs matrix vector multiplication for small matrices
# b gets overwritten
# multiplication is expressed as linear combination of columns of A
# optimized for column major format
# does that matter if all operands fit in cache?
# faster than Julia's built-in matvec for for square matrices
# of size 128 or less
  (m,n) = size(A)
  xm = length(x)
  bm = length(b)

  @assert n == xm
  @assert m == bm

  # overwrite b, first column only
  @inbounds begin
    @simd for i=1:m
      b[i] = x[1]*A[i, 1]
    end

    for i=2:n  # loop across columns
      @simd for j=1:m  # loop down columns
        b[j] += A[j,i]*x[i]
      end
    end

  end  # end begin inbounds

  return b

end   # end of smallmatvec! function


function smallmatvec{T, T2}(A::AbstractArray{T,2}, x::AbstractArray{T2, 1})
  (m,n) = size(A)
  T3 = promote_type(T, T2)
  b = Array(T3, m)
  smallmatvec!(A, x, b)
end


# do A.'*x = b
function smallmatTvec!(A::AbstractMatrix, x::AbstractVector, b::AbstractVector)
  (m,n) = size(A)
  xm = length(x)
  bm = length(b)
  
  @assert m == xm
  @assert n == bm

  @inbounds begin
    for i=1:n  # loop over columns of A
      # perform action for each column of A
      # over write each entry
      b[i] = A[1, i]*x[1]

      # accumulate over rest of column
      @simd for j=2:m
        b[i] += A[j, i]*x[j]
      end
    end

  end

  return b
end

function smallmatTvec{T, T2}(A::AbstractArray{T, 2}, x::AbstractArray{T2, 1})

  (m,n) = size(A)
  T3 = promote_type(T, T2)
  b = Array(T3, n)
  smallmatTvec!(A, x, b)
end


function smallmatmat!{T, T2, T3}(A::AbstractArray{T, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::AbstractArray{T3, 2})
# TODO: this comment needs to be a @doc
# multiplies matrix A by matrix x, writing the solution to matrix b
# both dimensions of A and the final dimension of x are used for looping
# the array sizes are not checked explicitly
# this uses the same technique as smallmatvec!, simply multiplying A by the columns
# of x repeatedly, without making any copies
# Faster than Julia's built-in mat-mat for matrices without m=n~=28

  (m,n) = size(A)
  (xn, p) = size(x)
  (bm, bn) = size(b)

  @assert n == xn
  @assert m == bm
  @assert p == bn

  @inbounds begin
    for k=1:p  # loop over the column vectors of x
      # overwrite b, first column only
      @simd for i=1:m
        b[i, k] = x[1, k]*A[i, 1]
      end

      for i=2:n  # loop across columns
        @simd for j=1:m  # loop down columns
          b[j, k] += A[j,i]*x[i, k]
        end
      end
    end

  end   # end begin inbounds

  return b

end     # end of smallmatmat! function


function smallmatmat{T, T2}(A::AbstractArray{T,2}, x::AbstractArray{T2, 2})
  (m,n) = size(A)
  (xn, p) = size(x)
  T3 = promote_type(T, T2)
  b = Array(T3, m, p)
  smallmatmat!(A, x, b)
end


function smallmatmatT!{T, T2, T3}(A::AbstractArray{T, 2},
                                  x::AbstractArray{T2, 2},
                                  b::AbstractArray{T3, 2})
# TODO: this comment needs to be a @doc
# multiplies A by x.', storing result in b
# the threshold (compared against a copying transpose of A*(x.') = b
# is m=n~=18 with bounds checking on
# without bounds checking, m=n~=24 appears to be the limit

  (p, xn) = size(x)
  (m,n) = size(A)
  (bm, bn) = size(b)

  @assert n == xn
  @assert m == bm
  @assert p == bn

  @inbounds begin
    # overwrite b
    for i=1:m
      a_i = A[i, 1]
      @simd for j=1:p
        b[i,j] = a_i*x[j,1]
      end
    end

    # add to b
    for j=2:n  # loop across remaining columns of A
      for i=1:m  # loop down a column
        a_i = A[i, j]
        # multiply this entry by row i of x.'
        @simd for k=1:p
          b[i,k] += a_i*x[k, j]
        end
      end
    end      

  end   # end begin inbounds

  return b

end     # end of smallmatmatT! function


function smallmatmatT{T, T2}(A::AbstractArray{T,2}, x::AbstractArray{T2, 2})
  (m,n) = size(A)
  (p, xn) = size(x)
  T3 = promote_type(T, T2)
  b = Array(T3, n, p)
  smallmatmatT!(A, x, b)
end


function checkZeroRows{T <: Number}(A::AbstractArray{T,2}, tol::AbstractFloat)
# checks each row of a matrix for zeros 
# 2d matrices only
# returns the integer number of zero rows, and a Bool
# array telling which rows have all zeros
  (m,n) = size(A)

  zero_mask = zeros(Bool, m)  # record whether each row has all zeros
  zero_row = zeros(Bool, n)  # record results for a row
  for i=1:m
    fill!(zero_row, false)
    for j=1:n
      if abs(A[i,j]) < tol  # if zero entry
        zero_row[j] = true
      end
    end  # end loop over column

    rowsum = sum(zero_row)
    zero_mask[i] = (rowsum == n)  # true if all zeros in row
  end  # end loop over rows

  return sum(zero_mask), zero_mask

end   # end of checkZeroRows function


function checkZeroColumns{T <: Number}(A::AbstractArray{T,2}, tol::AbstractFloat)
# TODO: @doc this
# checks each column of a matrix for zeros
# 2d matrices only
# returns integer number of zero rows, and a Bool
# array telling which rows have all zeros
  (m,n) = size(A)

  zero_mask = zeros(Bool, n)  # record whether each column has all zeros
  zero_row = zeros(Bool, m)  # record results for a column
  for i=1:n
    fill!(zero_row, false)
    for j=1:m
      if abs(A[j, i]) < tol  # if zero entry
        zero_row[j] = true
      end
    end  # end loop over row

    rowsum = sum(zero_row)
    zero_mask[i] = (rowsum == n)  # true if all zeros in row
  end  # end loop over columns

  return sum(zero_mask), zero_mask

end   # end checkZeroColumns function


function checkIdenticalColumns{T <: Number}(A::AbstractArray{T,2}, 
                                            colnum::Integer, 
                                            tol::AbstractFloat)
# TODO: @doc this
# checks which columns are identical to column number col
# returns number of column identical and an array of bools telling which ones
# does not count column colnum as being the same as itself

  (m,n) = size(A)
  cnt = 0
  
  col = view(A, :, colnum)
  is_same = zeros(Bool, n)
  
  for i=1:n  # loop over columns
  
    if i == colnum  # skip the specified column
      continue
    end
    col_i = view(A, :, i)
    diff_norm = norm(col - col_i)/m
  
    if diff_norm < tol
      is_same[i] = true
      cnt += 1
    end
  end     # end of for loop over columns
  
  return cnt, is_same

end     # end checkIdenticalColumns


function findLarge{T <: Number}(A::AbstractArray{T,2}, tol::AbstractFloat)
# TODO: @doc this

  (m,n) = size(A)
  cnt = 0
  for i=1:m
    for j=1:n
      if A[i,j] > tol
        println(i, ", ", j)
        cnt += 1
      end
    end
  end

  return cnt

end     # end findLarge


function checkSparseColumns{T <: Number, T2 <: Integer}(A::AbstractArray{T,2}, 
                            sparsity_bnds::AbstractArray{T2, 2}, 
                            tol::AbstractFloat)
# TODO: @doc this
# checks that all entries outside the range specified by sparsity_bnds
# are zero
# returns the number of columns with out of bounds entries and an array of
# bools specifying which ones

  (m,n) = size(A)
  out_of_bounds = zeros(Bool, n)
  
  for i=1:n  # loop over columns
    min = sparsity_bnds[1, i]
    max = sparsity_bnds[2, i]
  
    for j=1:(min -1)
      entry_j = A[j, i]
      if abs(entry_j) > tol
        out_of_bounds[i] = true
        println("entry ", j, ", ", i, " is non zero")
        break
      end
    end
  
    for j=(max+1):m
      entry_j = A[j, i]
      if abs(entry_j) > tol
        out_of_bounds[i] = true
        println("entry ", j, ", ", i, " is non zero")
        break
      end
    end
  
  end # end loop over columns
  
  cnt = sum(out_of_bounds)
  
  return cnt, out_of_bounds

end     # end of checkSparseColumns function


function checkSparseRows{T <: Number, T2 <: Integer}(A::AbstractArray{T,2}, 
                         sparsity_bnds::AbstractArray{T2, 2},  
                         tol::AbstractFloat)
# checks that all entries outside the range specified by sparsity_bnds
# are zero
# returns the number of rows with out of bounds entries and an array of
# bools specifying which ones

  (m,n) = size(A)
  out_of_bounds = zeros(Bool, m)
  
  for i=1:m  # loop over columns
    min = sparsity_bnds[1, i]
    max = sparsity_bnds[2, i]
  
    for j=1:(min -1)
      entry_j = A[i, j]
      if abs(entry_j) > tol
        out_of_bounds[i] = true
        println("entry ", i, ", ", j, " is non zero")
        break
      end
    end
  
    for j=(max+1):n
      entry_j = A[i, j]
      if abs(entry_j) > tol
        out_of_bounds[i] = true
        println("entry ", i, ", ", j, " is non zero")
        break
      end
    end
  
  end     # end loop over columns
  
  cnt = sum(out_of_bounds)
  
  return cnt, out_of_bounds

end     # end of checkSparseRows function


@doc """
### ODLCommonTools.isSymmetric

  This function checks if an array is symmetric or not, using
  the specified tolerance for comparing if two entries are
  equal.

  Inputs:
    A : an array to check for symmetry, must be possible to access all
        entries
    tol: tolerance for floating point equality

  Outputs:
    val: a Bool indicating if A is symmetric

"""->
function isSymmetric(A::AbstractArray, tol=1e-14)

  (m,n) = size(A)
  val = true
  for j=1:n
    for i=1:j
      val = val && (abs(A[i,j] - A[j,i]) < tol)
    end
  end

  return val
end

@doc """
ODLCommonTools.make_symmetric!

  This function copies the lower triangle of a matrix to the upper triangle,
  ensuring it is exactly symmetric.  It does not checking to determine if
  the matrix is close so symmetric beforehand

  Inputs/Outputs:
    A : a matrix of any type

"""->
function make_symmetric!(A::AbstractMatrix)
# make the matrix symmetrix, no questions asked
# the lower triangle of the matrix is copied to the upper

  for i=1:size(A,1)
    for j=1:(i-1)
      A[j, i] = A[i,j]
    end
  end

  return nothing
end

#----------------------------------------------------------
export FIFOQueue, front
import Base.push!, Base.pop!, Base.length, Base.isempty, 
       Base.resize!, Base.empty!
type FIFOQueue{T}
  s::Array{T, 1}  # array of values
  pos::Int  # current head of que (points to most recently inserted elements)
  fac::Float64  # factor by which to expand que when it runs out of space

  function FIFOQueue(; size_hint=0, fac=1.4)
    # size_hint = initial size of array
    arr = Array(T, size_hint)
    pos = 0
   return new(arr, pos, fac)
  end

end

function push!{T}(que::FIFOQueue{T}, val::T)
  que.pos += 1
  # check size
  println("que.pos = ", que.pos)
  println("que length = ", length(que.s))
  len = length(que.s)
  if que.pos > len
    println("resizing que")
    new_size = convert(Int, max(que.pos, floor(len*que.fac)))
    println("new_size = ", new_size)
    resize!(que.s, new_size)
  end

  println("que.pos = ", que.pos)
  # do insertion
  que.s[que.pos] = val

  return nothing
end

function pop!{T}(que::FIFOQueue{T})
  pos = que.pos
  val = que.s[pos]
  que.pos -= 1

  return val
end

function pop!{T}(que::FIFOQueue{T}, vals::AbstractArray{T, 1})
# remove last n elements, where n = length(vals)

  n = length(vals)
  for i=1:n
    vals[i] = pop!(que)
  end

  return nothing
end

# function front retrieve element without removing
function front{T}(que::FIFOQueue{T})
  return que.s[que.pos]
end

# function popN!  # remove multiple elements

function length(que::FIFOQueue)
  return que.pos
end

# for general collection compatability
endof(que::FIFOQueue) = length(que)

function isempty(que::FIFOQueue)
  return que.pos == 0
end

function empty!(que::FIFOQueue)
  n = length(que)

  for i=1:n
    pop!(que)
  end

  return nothing
end

function resize!(que::FIFOQueue, new_size)
  # change size of que, if doing so will not remove elements
  if new_size > que.pos
    resize!(que.s, new_size)
  else
    println(STDERR, "Warning: not resizing que")
  end

  return nothing
end


