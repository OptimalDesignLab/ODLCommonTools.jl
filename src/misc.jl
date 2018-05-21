# Description: miscellaneous functions

using ArrayViews

export rmfile, printbacktrace, smallmatvec!, smallmatvec, smallmatTvec!, 
        smallmatTvec, smallmatmat!, 
        smallmatmat, smallmatmatT!, smallmatmatT, smallmatTmat!, smallmatTmat,
        smallmatvec_revv!,
        checkZeroRows, 
        checkZeroColumns, checkIdenticalColumns, checkSparseColumns,
        checkSparseRows, findLarge, isSymmetric, make_symmetric!,
        getBranchName, getTimeString, isFieldDefined, get_parallel_fname,
        joinpath_ascii

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

function printbacktrace(file=STDOUT)
# print the current call stack
  bt = backtrace()
  s = sprint(io->Base.show_backtrace(io, bt))
  println(file,"backtrace: ", s)

  return s
end

#------------------------------------------------------------------------------
# Ax = b

"""
  Computes the product Ax = b

  This is a wrapper around different matrix-vector multiplication routines.
  It calls a special method for very small matrices and vectors, falling back
  to BLAS otherwise.  Performance tests indicate the check for which function
  to call is not that expensive

  **Inputs**
   
   * A: the matrix
   * x: the vector

  **Inputs/Outputs**

   * b: the output vector

  Aliasing restrictions: same as BLAS, b cannot alias anything
"""
function smallmatvec!{T, T2, T3}(A::AbstractArray{T,2}, 
           x::AbstractArray{T2,1}, b::AbstractArray{T3, 1})

  m, n = size(A)

  # extrapolating from the square matrix case, where my algorithm
  # is faster for m = n <= 8
  # for mixed float-int operations, the arrays need to be
  # reinterpratable
  if m*n <= 65
    smallmatvec_kernel!(A, x, b)
  else
    A_mul_B!(b, A, x)
#    LinAlg.BLAS.gemv!('N', 1.0, A, x, 0.0, b)
  end

  return b
end

# if either of the arguments that get reinerpreted are a ContiguousView,
# don't call BLAS
function smallmatvec!(A::ContiguousView{Complex128,2}, 
           x::AbstractArray{Float64,1}, b::ContiguousView{Complex128, 1})
# this is the case where we can't reinterpret as required by the BLAS interface

  smallmatvec_kernel!(A, x, b)

  return b
end


function smallmatvec!(A::ContiguousView{Complex128,2}, 
           x::AbstractArray{Float64,1}, b::AbstractArray{Complex128, 1})
# this is the case where we can't reinterpret as required by the BLAS interface

  smallmatvec_kernel!(A, x, b)

  return b
end

function smallmatvec!(A::AbstractArray{Complex128,2}, 
           x::AbstractArray{Float64,1}, b::ContiguousView{Complex128, 1})
# this is the case where we can't reinterpret as required by the BLAS interface

  smallmatvec_kernel!(A, x, b)

  return b
end

function smallmatvec!{P<:Array, P2<:Array}(A::ROContiguousView{Complex128,2, P}, 
           x::AbstractArray{Float64,1}, b::ContiguousView{Complex128, 1, P2})
# this is the case where we can't reinterpret as required by the BLAS interface

  smallmatvec_kernel!(A, x, b)

  return b
end


function smallmatvec!{P<:Array}(A::ROContiguousView{Complex128,2, P}, 
           x::AbstractArray{Float64,1}, b::AbstractArray{Complex128, 1})
# this is the case where we can't reinterpret as required by the BLAS interface

  smallmatvec_kernel!(A, x, b)

  return b
end





function smallmatvec_kernel!{T, T2, T3}(A::AbstractArray{T,2}, 
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

# reverse mode to back propigate b to x
function smallmatvec_revv!{T, T2, T3}(A::AbstractArray{T,2}, 
           x_bar::AbstractArray{T2,1}, 
           b_bar::AbstractArray{T3, 1})

  (m,n) = size(A)
  xm = length(x_bar)
  bm = length(b_bar)

  @assert n == xm
  @assert m == bm

  @inbounds begin
    # reverse mode always updates the output, so no need to special case
    # the i = 1 to overwrite it
    for i=1:n
      @simd for j=1:m
        x_bar[i] += b_bar[j]*A[j, i]
      end
    end

  end  # end begin inbounds

  return x_bar
end


function smallmatvec{T, T2}(A::AbstractArray{T,2}, x::AbstractArray{T2, 1})
  (m,n) = size(A)
  T3 = promote_type(T, T2)
  b = Array(T3, m)
  smallmatvec!(A, x, b)
end

#------------------------------------------------------------------------------
# A.'x = b


"""
  This function does A.'*x = b.  Note this uses the simple transpose of A, not
  the Hermetian transpose.

  This function call BLAS for large matrices and uses a special algorithm for
  small matrices.

  TODO: double check that the cutoff is correct.

  **Inputs**

   * A: the matrix
   * x: the vector

  **Inputs/Outputs**

   * b: vector overwritten with the result
"""
function smallmatTvec!(A::AbstractMatrix, x::AbstractVector, b::AbstractVector)
# currently this doesn't call BLAS for the complex-real case, so its ok to
# pass a ContiguousView to At_mul_B
  m, n = size(A)

  if m*n <= 65
    smallmatTvec_kernel!(A, x, b)
  else
    At_mul_B!(b, A, x)
  end

  return b
end


# do A.'*x = b
function smallmatTvec_kernel!(A::AbstractMatrix, x::AbstractVector, b::AbstractVector)
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

#------------------------------------------------------------------------------
# matrix-matrix multiply

"""
  Performs a matrix-matrix multiply.  Uses a special algorithm for small
  matrices, otherwise calls BLAS (when possible).

  **Inputs**

   * A: first matrix
   * x: second matrix

  **Inputs/Inputs**

   * x: output matrix, overwritten

   Aliasing restrictions: no aliasing allowed
"""
function smallmatmat!{T, T2, T3}(A::AbstractArray{T, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::AbstractArray{T3, 2})
  m, n = size(A)
  if m < 6 && n < 9
    smallmatmat_kernel!(A, x, b)
  else
    A_mul_B!(b, A, x)
  end

  return b
end

# don't call BLAS when reinterpret isn't supported

function smallmatmat!{T2}(A::ContiguousView{Complex128, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::ContiguousView{Complex128, 2})

  smallmatmat_kernel!(A, x, b)

  return b
end
 
function smallmatmat!{T2}(A::AbstractArray{Complex128, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::ContiguousView{Complex128, 2})

  smallmatmat_kernel!(A, x, b)

  return b
end
 
function smallmatmat!{T2}(A::ContiguousView{Complex128, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::AbstractArray{Complex128, 2})

  smallmatmat_kernel!(A, x, b)

  return b
end
 
function smallmatmat!{T2, P<:Array}(A::ROContiguousView{Complex128, 2, P}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::ContiguousView{Complex128, 2})

  smallmatmat_kernel!(A, x, b)

  return b
end
 
function smallmatmat!{T2, P<:Array}(A::ROContiguousView{Complex128, 2, P}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::AbstractArray{Complex128, 2})

  smallmatmat_kernel!(A, x, b)

  return b
end
 
function smallmatmat_kernel!{T, T2, T3}(A::AbstractArray{T, 2}, 
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

#------------------------------------------------------------------------------
# matrix-matrix transposed multiplication


"""
  Performs A*(x.') = b.  Uses a special algorithm for small matrices,
  otherwise calls BLAS

  **Inputs**

   * A: the first matrix
   * x: the matrix to be transposed

  **Inputs/Outputs**

   * b: output matrix, overwritten
"""
function smallmatmatT!{T, T2, T3}(A::AbstractArray{T, 2},
                                  x::AbstractArray{T2, 2},
                                  b::AbstractArray{T3, 2})
  m, n = size(A)
  if m < 7 && n < 10
    smallmatmatT_kernel!(A, x, b)
  else
    A_mul_Bt!(b, A, x)
  end

  return b
end

# don't call BLAS when reinterpret isn't supported

function smallmatmatT!{T2}(A::ContiguousView{Complex128, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::ContiguousView{Complex128, 2})

  smallmatmatT_kernel!(A, x, b)

  return b
end
 
function smallmatmatT!{T2}(A::AbstractArray{Complex128, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::ContiguousView{Complex128, 2})

  smallmatmatT_kernel!(A, x, b)

  return b
end
 
function smallmatmatT!{T2}(A::ContiguousView{Complex128, 2}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::AbstractArray{Complex128, 2})

  smallmatmatT_kernel!(A, x, b)

  return b
end
 
function smallmatmatT!{T2, P<:Array}(A::ROContiguousView{Complex128, 2, P}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::ContiguousView{Complex128, 2})

  smallmatmatT_kernel!(A, x, b)

  return b
end
 
function smallmatmatT!{T2, P<:Array}(A::ROContiguousView{Complex128, 2, P}, 
                                 x::AbstractArray{T2, 2}, 
                                 b::AbstractArray{Complex128, 2})

  smallmatmatT_kernel!(A, x, b)

  return b
end
 
function smallmatmatT_kernel!{T, T2, T3}(A::AbstractArray{T, 2},
                                  x::AbstractArray{T2, 2},
                                  b::AbstractArray{T3, 2})
# TODO: this comment needs to be a @doc
# multiplies A by x.', storing result in b

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

#------------------------------------------------------------------------------
# matrix transposed-matrix multiply

"""
  Performs A.'*x = b.  Uses a special algorithm for small matrices, otherwise
  calls BLAS

  **Inputs**

   * A: the matrix to be transposed
   * x: the second matrix

  **Inputs/Outputs**

   * b: output matrix, overwritten

  Aliasing restrictions: no aliasing
"""
function smallmatTmat!{T, T2, T3}(A::AbstractMatrix{T},
                                         x::AbstractMatrix{T2},
                                         b::AbstractMatrix{T3})
  m, n = size(A)
  if m < 10 && n < 7
    smallmatTmat_kernel!(A, x, b)
  else
    At_mul_B!(b, A, x)
  end

  return b
end

# this doesn't call BLAS in the Complex128-Float64 case, so no need to deal
# with reinterpret issues

#TODO: performance test
function smallmatTmat_kernel!{T, T2, T3}(A::AbstractMatrix{T}, x::AbstractMatrix{T2},
                                  b::AbstractMatrix{T3})

  n, m = size(A)
  xn, p = size(x)
  bm, bn = size(b)

  @assert n == xn
  @assert m == bm
  @assert p == bn

  # overwrite b
  fill!(b, 0)

  @inbounds for i=1:p  # columns of x
    for j=1:m  # columns of A
      # zero b here?
      @simd for k=1:n
        b[j, i] += A[k, j]*x[k, i]
      end
    end
  end

  return nothing
end

function smallmatTmat{T, T2}(A::AbstractMatrix{T}, x::AbstractMatrix{T2})

  n, m = size(A)
  n, p = size(x)
  T3 = promote_type(T, T2)
  b = zeros(m, p)
  smallmatTmat!(A, x, b)
  return b
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
  tail::Int  # current tail of que (points to most recently inserted elements)
  head::Int  # points to least recently inserted elements (next to pop)
  fac::Float64  # factor by which to expand que when it runs out of space

  function FIFOQueue(; size_hint=1, fac=1.4)
    # size_hint = initial size of array
    size_hint = max(size_hint, 1)
    arr = Array(T, size_hint)
    tail = 0
    head = 1
   return new(arr, tail, head, fac)
  end

end

function push!{T}(que::FIFOQueue{T}, val::T)
  # check size
  len = length(que.s)
  if (que.tail + 1) > len
    resize_que!(que)
  end

  que.tail += 1
  
  # do insertion
  que.s[que.tail] = val

  return nothing
end

function pop!{T}(que::FIFOQueue{T})
  # add check that head >= tail ?
  pos = que.head
  val = que.s[pos]
  que.head += 1

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
  return que.s[que.head]
end

# function popN!  # remove multiple elements

function length(que::FIFOQueue)
  return que.tail - que.head + 1
end

# for general collection compatability
endof(que::FIFOQueue) = length(que)

function isempty(que::FIFOQueue)
  return que.head > que.tail
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
  if new_size > que.tail
    resize!(que.s, new_size)
  else
    println(STDERR, "Warning: not resizing que")
  end

  return nothing
end

# for resizing the que when it runs out of space
function resize_que!(que::FIFOQueue)

  # decide if shifting the queue or resizing it is better
  if que.head > cld(length(que), 10)
    pos = 1
    for i=que.head:que.tail
      que.s[pos] = que.s[i]
      pos += 1
    end
    que.head = 1
    que.tail = pos - 1
  else  # resize
    new_size = convert(Int, ceil(length(que)*que.fac))
    resize!(que.s, new_size)
  end
end


function getBranchName(dir=pwd())
# get the name of the current branch of the git repo in the specified directory
  nm = readall(`git rev-parse --abbrev-ref HEAD`)
  return nm[1:end-1]  # remove newline

return end

function getTimeString()
  t = now()
  y = Dates.year(t)
  m = Dates.month(t)
  d = Dates.day(t)
  h = Dates.hour(t)
  minutes = Dates.minute(t)

  return "$y-$m-$d $h:$minutes"
end
@doc """
  
  Wrapper and generalization of Base.isdefined.
  This function takes an object and one or more symbols and checks if
  the fields of the object with those names are defined, but throws an 
  exception if the symbol is not a field of the type.  Returns true 
  iff all field names are defined.
"""
function isFieldDefined(obj, req_fieldnames...)

  if length(req_fieldnames) == 0
    throw(ErrorException("must check at least one fieldname"))
  end

  
  fnames = fieldnames(obj)
  obj_type = typeof(obj)
  alldefined = true
  for i=1:length(req_fieldnames)
    fname_i = symbol(req_fieldnames[i])  # convert to symbol if possible

    if !(fname_i in fnames)
      throw(ErrorException("fieldname $fname_i is not a field of $obj_type"))
    end

    alldefined = alldefined && isdefined(obj, fname_i)
  end

  return alldefined
end


"""
  Given a file name (including extension), adds an _comm_rank, where comm_rank
  is an MPI communicator rank, to the end of the name, before the extension.
  For example, foo.dat -> foo_0.dat

  Also works if the file name does not have an extension:

  foo -> foo_0.dat

  Inputs:
    fname: original file name, including extension
    comm_rank: communicator rank

  Outputs:
    fname_stub: new file name, including extension
"""
function get_parallel_fname(fname::ASCIIString, comm_rank)

  # figure out where file extension starts
  len = length(fname)
  sep_loc = 0  # location of separator
  for i=len:-1:1
    if fname[i] == '.'
      sep_loc = i
    end
    if fname[i] == '/'  # linux only
      break
    end
  end

  if sep_loc == 0  # fname does not have an extension
    fname_stub = string(fname, "_", comm_rank)
  else
    fname_stub = fname[1:(sep_loc-1)]  # grab the name before the extension
    fname_stub *= "_$comm_rank"  # add the comm_rank
    fname_stub *= fname[sep_loc:end]
  end

  return fname_stub
end

"""
  Returns the file name and extension as separate strings.
  Ex. output.dat -> output, .dat

  If the file does not have an extension the second string will be empty
  **Input**

   * fname: file name to split

  **Outputs**

   * fstem: file name, without extension
   * fext: file extension (including the period)
"""
function split_fname(fname::AbstractString)

  sep_loc = 0
  for i=length(fname):-1:1
    if fname[i] == '.'
      sep_loc = i
    end
  end

  if sep_loc == 0  # no extension
    fstem = copy(fname)
    fext = ""
  else
    fstem = fname[1:(sep_loc-1)]
    fext = fname[sep_loc:end]
  end


  return fstem, fext
end

"""
  Wrapper around joinpath() that always returns an ASCIIString.
"""
function joinpath_ascii(str::ASCIIString...)

  return ASCIIString(joinpath(str...))
end

"""
  Prepends the given string to another string as though it is a Linux path
  (ie. using a colon delimiter)

  Inputs:

    path: the existing path, can be an empty string
    new_entry: string to add to the path (should not contain the delimiter)

  Outputs:

    new_path: the results of prepend operation
"""
function prepend_path(path::AbstractString, new_entry::AbstractString)

  if path == ""
    new_path = new_entry
  else
    new_path = string(new_entry, ":", path)
  end

  return new_path
end

"""
  Like [`prepend_path`](@ref), but adds the new entry to the end rather than
  the beginning
"""
function append_path(path::AbstractString, new_entry::AbstractString)

  if path == ""
    new_path = new_entry
  else
    new_path = string(path, ":", new_entry)
  end

  return new_path
end



