# a read-only wrapper for an AbstractArray
"""
  A read only wrapper type for AbstractArrays. getindex is defined for
  this function but setindex! is not.
"""
immutable ROView{T, N, P <: AbstractArray} <: AbstractArray{T, N}
  parent::P
end

"""
  Outer constructor
"""
function ROView{T, N}(A::AbstractArray{T, N})
  return ROView{T, N, typeof(A)}(A)
end


# define some functions by delegating to parent array
import Base: size, length, getindex, parent, convert, linearindexing

size(A::ROView) = size(A.parent)
length(A::ROView) = length(A.parent)

# it appears there is no splatting penalty in this case
getindex(A::ROView, idx...) = getindex(A.parent, idx...)

parent(A::ROView) = A.parent

linearindexing(A::ROView) = linearindexing(A.parent)

"""
  Create a sview and then make a read only view of it
"""
@inline ro_sview(A::Array, idx...) = ROView(sview(A, idx...))
@inline ro_sview(A::ROView, idx...) = ROView(sview(parent(A), idx...))
@inline ro_sview(A::ArrayView, idx...) = ROView(sview(A, idx...))

# useful aliases
typealias ROVector{T} ROView{T, 1}
typealias ROMatrix{T} ROView{T, 2}
typealias ROArray{T, N} ROView{T, N}

convert{T, T2}(::Type{Array{T}}, x::ROView{T2}) = ROView(convert(Array{T}, parent(x)))

copy!{T, N, T2}(dest::AbstractArray{T, N}, src::ROView{T2, N}) = copy!(dest, src.parent)
