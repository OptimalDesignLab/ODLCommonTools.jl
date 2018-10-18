# a read-only wrapper for an AbstractArray
# because julia's type system is not sufficiently expressive, this will
# only work for arrays (needs to be subtype of DenseArray in order 
# to call BLAS functions

import ArrayViews: aview, unsafe_aview, Subs

ContiguousArrays{T, N} =   Union{Array{T, N}, ContiguousView{T, N}, UnsafeContiguousView{T, N}}

const AllArrayViews = Union{ArrayView, UnsafeArrayView}

"""
  A read only wrapper type for AbstractArrays. getindex is defined for
  this function but setindex! is not.
"""
struct ROView{T, N, P <: AbstractArray} <: DenseArray{T, N}
  parent::P

  function ROView{T, N, P}(a::ContiguousArrays) where {T, N, P}
    return new(a)
  end
end

"""
  Outer constructor
"""
function ROView(A::ContiguousArrays{T, N}) where {T, N}
  return ROView{T, N, typeof(A)}(A)
end


# define some functions by delegating to parent array
import Base: size, length, getindex, parent, convert, IndexStyle, convert, unsafe_convert, pointer

size(A::ROView) = size(A.parent)
length(A::ROView) = length(A.parent)

# it appears there is no splatting penalty in this case
@inline getindex(A::ROView, idx...) = getindex(A.parent, idx...)

parent(A::ROView) = A.parent

IndexStyle(A::ROView) = IndexStyle(A.parent)

pointer(A::ROView) = pointer(A.parent)

unsafe_convert(::Type{Ptr{T}}, A::ROView) where {T} = pointer(A)

"""
  Create a sview and then make a read only view of it
"""
@inline ro_sview(A::Array, idx...) = ROView(sview(A, idx...))
@inline ro_sview(A::ROView, idx...) = ROView(sview(parent(A), idx...))
@inline ro_sview(A::AllArrayViews, idx...) = ROView(sview(A, idx...))

# enable taking a sview of a ro_sview
# this must return an ro_sview to preserve the read-only property
# up to 5 argument constructors have to be listed explicitly to avoid ambiguity
# problem with Arrayviews (a varargs is less specific than a positional argument
# and ROView us a subtype of DenseArray

@inline aview(A::ROView, i1::Subs) = ROView(aview(parent(A), i1))
@inline aview(A::ROView, i1::Subs, i2::Subs) = ROView(aview(parent(A), i1, i2))
@inline aview(A::ROView, i1::Subs, i2::Subs, i3::Subs) = ROView(aview(parent(A), i1, i2, i3))
@inline aview(A::ROView, i1::Subs, i2::Subs, i3::Subs, i4::Subs) = ROView(aview(parent(A), i1, i2, i3, i4))
@inline aview(A::ROView, i1::Subs, i2::Subs, i3::Subs, i4::Subs, i5::Subs) = ROView(aview(parent(A), i1, i2, i3, i4, i5))
@inline aview(A::ROView, i1::Subs, i2::Subs, i3::Subs, i4::Subs, i5::Subs, I::Subs...) = ROView(aview(parent(A), i1, i2, i3, i4, i5, I...))


@inline unsafe_aview(A::ROView, i1::Subs) = ROView(unsafe_aview(parent(A), i1))
@inline unsafe_aview(A::ROView, i1::Subs, i2::Subs) = ROView(unsafe_aview(parent(A), i1, i2))
@inline unsafe_aview(A::ROView, i1::Subs, i2::Subs, i3::Subs) = ROView(unsafe_aview(parent(A), i1, i2, i3))
@inline unsafe_aview(A::ROView, i1::Subs, i2::Subs, i3::Subs, i4::Subs) = ROView(unsafe_aview(parent(A), i1, i2, i3, i4))
@inline unsafe_aview(A::ROView, i1::Subs, i2::Subs, i3::Subs, i4::Subs, i5::Subs) = ROView(unsafe_aview(parent(A), i1, i2, i3, i4, i5))
@inline unsafe_aview(A::ROView, i1::Subs, i2::Subs, i3::Subs, i4::Subs, i5::Subs, I::Subs...) = ROView(unsafe_aview(parent(A), i1, i2, i3, i4, i5, I...))

# useful aliases
ROVector{T, P} =  ROView{T, 1, P}
ROMatrix{T, P} =  ROView{T, 2, P}
ROArray{T, N, P} =  ROView{T, N, P}
ROContiguousView{T, N, P} =  ROView{T, N, ContiguousView{T, N, P}}


convert(::Type{Array{T}}, x::ROView{T2}) where {T, T2} = ROView(convert(Array{T}, parent(x)))

copy!(dest::AbstractArray{T, N}, src::ROView{T2, N}) where {T, N, T2} = copy!(dest, src.parent)


import Base.reinterpret
"""
  Reinterpret an ArrayViews.UnsafeContiguousView (needed by BLAS wrappers).

  It is not possible, in general, to reinterpert a regular ContiguousView, so
  the small linear algebra functions will have to avoid calling BLAS routines in
  that case
"""
function reinterpret(::Type{T}, a::UnsafeContiguousView{S},
                     dims::NTuple{N, Int}) where {T, S, N}


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
  new_arr = UnsafeContiguousView(pnew, dims)

  return new_arr
end

function reinterpret(::Type{T}, a::UnsafeContiguousView{S}) where {T, S}
# Blas mat-vec wrapper uses this

  size_in = sizeof(S)
  size_out = sizeof(T)

  # if this isn't exact, the length check later will error
  first_dim = size(a, 1)*div(size_in, size_out)
  new_dims = (first_dim, size(a)[2:end]...)

  reinterpret(T, a, new_dims)
end


function reinterpret(::Type{T}, a::UnsafeContiguousView{T}) where T

  return a
end

function reinterpret(::Type{T}, a::ContiguousView{T}) where T

  return a
end

function reinterpret(::Type{T}, a::ROView{S}, dims::NTuple{N, Int}) where {T, S, N}
# reinterpret the parent and wrap it in a ROView

  return ROView(reinterpret(T, parent(a), dims))
end

function reinterpret(::Type{T}, a::ROView{T}) where T

  return a
end
#=
function reinterpret{T}(::Type{T}, a::Union{ContiguousView{T}, UnsafeContiguousView{T}})
# it is always possible to reinterpret to the same type

  return a
end
=#

#TODO: reinterpret a ROView as well
