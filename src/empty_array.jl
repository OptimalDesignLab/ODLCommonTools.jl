# an (stack-allocatable) array for which setindex and getindex are no-ops
# This is useful for unneeded outputs

import Base: size, getindex, setindex!, IndexStyle

"""
  Special array type that errors if indexed.  
  This is useful for the idiom in C of passing a null pointer to a function
  to signal that an argument is was not supplied

  This array can be stack-allocated.  It supportes the standard Array
  constructors.
"""
struct EmptyArray{T, N} <: DenseArray{T, N}
  dims::NTuple{N, Int}
  offset::Int  # compatability with ArrayViews

  function EmptyArray{T, N}(dims::NTuple{N, Int}, offset::Integer) where {T, N}
    return new(dims, offset)
  end

  function EmptyArray{T, N}(dims::NTuple{N, Int}) where {T, N}
    return new(dims, 0)
  end
end

IndexStyle(::EmptyArray) = IndexLinear()
size(A::EmptyArray) = A.dims

function getindex(A::EmptyArray{T}, i::Int) where {T}
  return zero(T)  # this makes += operator work
end


function setindex!(A::EmptyArray, v, i::Int)
  return A  # I think this is what Base does for the return value
end

# construct from individual dimensions
# These constructors can be called as EmptyArray{T}(i, j) (which is not
# documented in Julia).
function EmptyArray{T}(i::Integer) where {T}
  return EmptyArray{T, 1}( (Int(i), ))
end

function EmptyArray{T}(i::Integer, j::Integer) where {T}
  return EmptyArray{T, 2}( (Int(i), Int(j)))
end

function EmptyArray{T}(i::Integer, j::Integer, k::Integer) where {T}
  return EmptyArray{T, 3}( (Int(i), Int(j), Int(k)) )
end

function EmptyArray{T}(i::Integer, j::Integer, k::Integer, l::Integer
                         ) where {T}
  return EmptyArray{T, 4}( (Int(i), Int(j), Int(k), Int(l)) )
end

function EmptyArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer) where {T}
  return EmptyArray{T, 5}( (Int(i), Int(j), Int(k), Int(l), Int(m)) )
end

function EmptyArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer, n::Integer) where {T}
  return EmptyArray{T, 6}( (Int(i), Int(j), Int(k), Int(l), Int(m), Int(n)) )
end

function EmptyArray{T}(i::Integer, j::Integer, k::Integer, l::Integer,
                          m::Integer, n::Integer, o::Integer) where {T}
  return EmptyArray{T, 7}( (Int(i), Int(j), Int(k), Int(l), Int(m), Int(n), Int(o)) )
end

# NTuple constructor
function EmptyArray{T}(d::NTuple{T, Int})
  return EmptyArray{T, N}(d)
end


#------------------------------------------------------------------------------
# compatability with sview (limited to contiguous views, 2D parent array)

import ODLCommonTools.sview

function sview(A::EmptyArray{T, 2}, rng::UnitRange, i::Integer) where {T}
  dims = (Int(length(rng)), )
  offset = rng[1] - 1
  return EmptyArray{T, 1}(dims, offset)
end

function sview(A::EmptyArray{T, 2}, idx::Colon, i::Integer) where {T}

  dims = (Int(A.dims[1]),)
  return EmptyArray{T, 1}(dims)
end


