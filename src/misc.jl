using ArrayViews

export rmfile, printbacktrace, smallmatvec!, checkZeroRows, checkZeroColumns, checkIdenticalColumns, checkSparseColumns

@doc """
 ### Tools rmfile

 Removes the file if it exists.  This is only necessary because it is an error
   to remove a file that does not exist.  Does not work on directories

   Arguments:
     fname::AbstractString  :  name of file

"""->
function rmfile(fname::AbstractString)

  if isfile(fname)
    rm(fname)
  end

  return nothing
end



function printbacktrace()
# print the current call stack
  bt = backtrace()
  s = sprint(io->Base.show_backtrace(io, bt))
  println("backtrace: ", s)

  return nothing
end


function smallmatvec!{T, T2, T3}(A::AbstractArray{T,2}, x::AbstractArray{T2,1}, b::AbstractArray{T3, 1})
# performs matrix vector multiplication for small matrices
# b gets overwritten
# multiplication is expressed as linear combination of columns of A
# optimized for column major format
# does that matter if all operands fit in cache?
(m,n) = size(A)

# overwrite b, first column only
for i=1:m
  b[i] = x[1]*A[i, 1]
end

for i=2:n  # loop across columns
  for j=1:m  # loop down columns
    b[j] += A[j,i]*x[i]
  end
end

  return nothing
end

function checkZeroRows{T <: Number}(A::AbstractArray{T,2}, tol::FloatingPoint)
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

end

function checkZeroColumns{T <: Number}(A::AbstractArray{T,2}, tol::FloatingPoint)
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

end

function checkIdenticalColumns{T <: Number}(A::AbstractArray{T,2}, colnum::Integer, tol::FloatingPoint)
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
end

return cnt, is_same

end

function findLarge{T <: Number}(A::AbstractArray{T,2}, tol::FloatingPoint)
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
end

function checkSparseColumns{T <: Number, T2 <: Integer}(A::AbstractArray{T,2}, sparsity_bnds::AbstractArray{T2, 2},  tol::FloatingPoint)
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

#  if fail_flag
#    continue  # skip to next column
#  end

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

end


function checkSparseRows{T <: Number, T2 <: Integer}(A::AbstractArray{T,2}, sparsity_bnds::AbstractArray{T2, 2},  tol::FloatingPoint)
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

#  if fail_flag
#    continue  # skip to next column
#  end

  for j=(max+1):n
    entry_j = A[i, j]
    if abs(entry_j) > tol
      out_of_bounds[i] = true
      println("entry ", i, ", ", j, " is non zero")
      break
    end
  end

end # end loop over columns

cnt = sum(out_of_bounds)

return cnt, out_of_bounds

end






export FIFOQueue, front
import Base.push!, Base.pop!, Base.length, Base.isempty, Base.resize!, Base.empty!
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






