# an (stack-allocatable) array that errors if you try to index it
# This is useful for the idiom in C of passing a null pointer to a function
# to signal that an argument is was not supplied

import Base: size, getindex, setindex!

"""
  Special array type that errors if indexed.  
  This is useful for the idiom in C of passing a null pointer to a function
  to signal that an argument is was not supplied

  This array can be stack-allocated.  It supportes the standard Array
  constructors.
"""
struct NullArray{T, N} <: DenseArray{T, N}
  dims::NTuple{N, Int}
  offset::Int  # compatability with ArrayViews

  function NullArray{T, N}(dims::NTuple{N, Int}, offset::Integer) where {T, N}
    return new(dims, offset)
  end

  function NullArray{T, N}(dims::NTuple{N, Int}) where {T, N}
    return new(dims, 0)
  end
end

size(A::NullArray) = A.dims

function getindex(A::NullArray{T, 1}, i::Int) where {T}
  throw(ErrorException("getindex! not defind for NullArray"))
end


function setindex!(A::NullArray, v, i::Int)
  throw(ErrorException("setindex! not defind for NullArray"))
end

# construct from individual dimensions
# These constructors can be called as NullArray{T}(i, j) (which is not
# documented in Julia).
function NullArray{T}(i::Integer) where {T}
  return NullArray{T, 1}( (Int(i), ))
end

function NullArray{T}(i::Integer, j::Integer) where {T}
  return NullArray{T, 2}( (Int(i), Int(j)))
end

function NullArray{T}(i::Integer, j::Integer, k::Integer) where {T}
  return NullArray{T, 3}( (Int(i), Int(j), Int(k)) )
end

function NullArray{T}(i::Integer, j::Integer, k::Integer, l::Integer
                         ) where {T}
  return NullArray{T, 4}( (Int(i), Int(j), Int(k), Int(l)) )
end

function NullArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer) where {T}
  return NullArray{T, 5}( (Int(i), Int(j), Int(k), Int(l), Int(m)) )
end

function NullArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer, n::Integer) where {T}
  return NullArray{T, 6}( (Int(i), Int(j), Int(k), Int(l), Int(m), Int(n)) )
end

function NullArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer, n::Integer, o::Integer) where {T}
  return NullArray{T, 7}( (Int(i), Int(j), Int(k), Int(l), Int(m), Int(n), Int(o)) )
end

# NTuple constructor
function NullArray{T}(d::NTuple{T, Int})
  return NullArray{T, N}(d)
end


#------------------------------------------------------------------------------
# compatability with sview (limited to contiguous views, 2D parent array)

import ODLCommonTools.sview

function sview(A::NullArray{T, 2}, rng::UnitRange, i::Integer) where {T}
  dims = (Int(length(rng)), )
  offset = rng[1] - 1
  return NullArray{T, 1}(dims, offset)
end

function sview(A::NullArray{T, 2}, idx::Colon, i::Integer) where {T}

  dims = (Int(A.dims[1]),)
  return NullArray{T, 1}(dims)
end


