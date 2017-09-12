# a read-only wrapper for an AbstractArray
# because julia's type system is not sufficiently expressive, this will
# only work for arrays (needs to be subtype of DenseArray in order 
# to call BLAS functions

typealias ROTypes{T, N}  Union{Array{T, N}, ContiguousView{T, N}, UnsafeContiguousView{T, N}}

"""
  A read only wrapper type for AbstractArrays. getindex is defined for
  this function but setindex! is not.
"""
immutable ROView{T, N, P <: AbstractArray} <: DenseArray{T, N}
  parent::P

  function ROView(a::ROTypes)
    return new(a)
  end
end

"""
  Outer constructor
"""
function ROView{T, N}(A::ROTypes{T, N})
  return ROView{T, N, typeof(A)}(A)
end


# define some functions by delegating to parent array
import Base: size, length, getindex, parent, convert, linearindexing, convert, unsafe_convert, pointer

size(A::ROView) = size(A.parent)
length(A::ROView) = length(A.parent)

# it appears there is no splatting penalty in this case
@inline getindex(A::ROView, idx...) = getindex(A.parent, idx...)

parent(A::ROView) = A.parent

linearindexing(A::ROView) = linearindexing(A.parent)

pointer(A::ROView) = pointer(A.parent)

unsafe_convert{T}(::Type{Ptr{T}}, A::ROView) = pointer(A)

"""
  Create a sview and then make a read only view of it
"""
@inline ro_sview(A::Array, idx...) = ROView(sview(A, idx...))
@inline ro_sview(A::ROView, idx...) = ROView(sview(parent(A), idx...))
@inline ro_sview(A::ArrayView, idx...) = ROView(sview(A, idx...))

# useful aliases
typealias ROVector{T, P} ROView{T, 1, P}
typealias ROMatrix{T, P} ROView{T, 2, P}
typealias ROArray{T, N, P} ROView{T, N, P}
typealias ROContiguousView{T, N, P} ROView{T, N, ContiguousView{T, N, P}}


convert{T, T2}(::Type{Array{T}}, x::ROView{T2}) = ROView(convert(Array{T}, parent(x)))

copy!{T, N, T2}(dest::AbstractArray{T, N}, src::ROView{T2, N}) = copy!(dest, src.parent)


import Base.reinterpret
"""
  Reinterpret an ArrayViews.UnsafeContiguousView (needed by BLAS wrappers).

  It is not possible, in general, to reinterpert a regular ContiguousView, so
  the small linear algebra functions will have to avoid calling BLAS routines in
  that case
"""
function reinterpret{T, S, N}(::Type{T}, a::UnsafeContiguousView{S},
                              dims::NTuple{N, Int})


# T is the output element type
# S is the input element type
# N is the output dimensionality

  # print the same warnings as the Base code
  # TODO: make sure these don't have a performance impact
  if !isbits(T)
    throw(ArgumentError("cannot reinterpret array to non bitstype $(T)"))
  end

  if !isbits(S)
    throw(ArgumentError("cannot reinterpret array from non bitstype $(S)"))
  end

  # ensure size in memory is same
  size_in = prod(size(a))*sizeof(S)
  len_out = prod(dims)
  size_out = len_out*sizeof(T)

  if size_in != size_out
    throw(DimensionMismatch("input and output sizes must be consistent"))
  end

  pnew = Ptr{T}(a.ptr)
  new_arr = UnsafeContiguousView{T, N}(pnew, len_out, dims)

  return new_arr
end

function reinterpret{T, S}(::Type{T}, a::UnsafeContiguousView{S})
# Blas mat-vec wrapper uses this

  size_in = sizeof(S)
  size_out = sizeof(T)

  # if this isn't exact, the length check later will error
  first_dim = size(a, 1)*div(size_in, size_out)
  new_dims = (first_dim, size(a)[2:end]...)

  reinterpret(T, a, new_dims)
end


function reinterpret{T}(::Type{T}, a::UnsafeContiguousView{T})

  return a
end

function reinterpret{T}(::Type{T}, a::ContiguousView{T})

  return a
end

function reinterpret{T, S, N}(::Type{T}, a::ROView{S}, dims::NTuple{N, Int})
# reinterpret the parent and wrap it in a ROView

  return ROView(reinterpret(T, parent(a), dims))
end

function reinterpret{T}(::Type{T}, a::ROView{T})

  return a
end
#=
function reinterpret{T}(::Type{T}, a::Union{ContiguousView{T}, UnsafeContiguousView{T}})
# it is always possible to reinterpret to the same type

  return a
end
=#

#TODO: reinterpret a ROView as well
