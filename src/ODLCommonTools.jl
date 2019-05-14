__precompile__(true)
# Name: ODLCommonTools
# Description: Common Tools for ODL
#   Performs the function of forward declaring abstract types

module ODLCommonTools
using ArrayViews
import ArrayViews.view  # use ArrayViews rather than Base.view
import Base.==

include("topo.jl")

import Base.show
import Base.isless
import Base.copy
import Base.copy!
export AbstractSolutionData, AbstractParamType, Abstract3DArray, Abstract4DArray
export AbstractMesh, AbstractCGMesh, AbstractDGMesh, AbstractOperator,
       AbstractSharedFaceData
export Boundary, Interface, getElementL, getFaceL, NullBoundary, NullInterface
export AbstractFunctional, AbstractIntegralFunctional,
       AbstractBoundaryFunctional, setupFunctional, getParallelData,
       getParallelDataString, getParallelDataEnum, CompositeFunctional,
       PARALLEL_DATA_NONE, PARALLEL_DATA_FACE, PARALLEL_DATA_ELEMENT
export Boundary
export Interface
export BCType, BCType_revm, BCType_revq, SRCType, FluxType, FluxType_revm,
       FluxType_diff, FluxType_revq, FunctionalType
export calcDiffElementArea
export functorThatErrors, functorThatErrors_revm

# topo.jl
export ElementTopology3, ElementTopology2, ElementTopology

# eqn_copy.jl
export copyForMultistage!

# eqn_deepcopy.jl
export eqn_deepcopy, eqn_deepcopy_fields

export ROView, ro_sview, ROVector, ROMatrix, ROArray

# misc.jl
export prepend_path, append_path, split_fname

# io.jl
export write_binary, read_binary!, writeSolutionFiles, readSolutionFiles

# lapack.jl
export getrf!, getrs2!, BlasInt, laswp!, applyIpiv!, solve_suitesparse,
       UMFPACK_A, UMFPACK_At, UMFPACK_Aat

export sview

# null_array.jl
export NullArray

# empty_array.jl
export EmptyArray

"""

  This abstract type is the supertype for all the objects that store the 
  solution data. Every physics module should implement its own subtype.

  **Static parameters**

   * Tsol: datatype of solution variables
   * Tres: datatype of the mesh variables

  See the [AbstractSolutionData](@ref) for the description of everything this
  type must implement.

"""
abstract type AbstractSolutionData{Tsol, Tres} end

"""
  This abstract type is the supertype for all mesh objects.  Every interface to
  a mesh software should define its own implementation.

  **Static parameters**

   * Tmsh: datatype of the mesh data (coordinates, mapping to/from parametric
          space, mapping jacobian).

  See the [AbstractMesh](@ref) for the description of everything this
  type must implement.

"""
abstract type AbstractMesh{Tmsh} end

@doc """
### ODLCommonTools.AbstractCGMesh

  The abstrac type is the supertype of all continuous Galerkin meshes
"""->
abstract type AbstractCGMesh{Tmsh} <: AbstractMesh{Tmsh} end

@doc """
### ODLCommonTools.AbstractDGGMesh

  The abstract type is the supertype of all discontinuous Galerkin meshes
"""->
abstract type AbstractDGMesh{Tmsh} <: AbstractMesh{Tmsh} end

"""
  Abstract type for all discretization operators (SBP, FE, etc.)
"""
abstract type AbstractOperator{T} end

"""
  This abstract type is the supertype for all Param objects, which hold values 
  needed for the computation in a place that is fast to access.

   The Param type is also useful for dispatching to low level functions which 
   the AbstractSolutionData might not be passed (depending on the organization 
   of the physics module.

  **Static Parameters**:
   
   * Tdim: the dimensionality of the equation being solved (2d or 3d usually)
"""
abstract type AbstractParamType{Tdim} end

@doc """
### ODLCommonTools.Abstract3DArray

  A typealias useful for specify a 3 dimensional AbstractArray without 
  specifying the element datatype.

"""->
const Abstract3DArray{T} = AbstractArray{T, 3}

"""
  A typealias for any 4D array.  Element is the static parameter
"""
const Abstract4DArray{T} =  AbstractArray{T, 4}

"""
Abstract supertype of all SharedFaceData implemenations, used to storing
the data needed for parallel communication

  **Static Parameters**

   * Tsol: the datatype of the solution variables
"""
abstract type AbstractSharedFaceData{Tsol} end



function show(io::IO, obj::AbstractMesh)
  println(io, "$(typeof(obj)) with $(obj.numEl) p = $(obj.order) elements")
end

function show(io::IO, obj::AbstractSolutionData)
  println(io, "$(typeof(obj)) with $(size(obj.q, 3)) elements")
end

function show(io::IO, obj::AbstractOperator)
  println(io, "degree $(obj.degree) $(typeof(obj))")
end


"""

  Copies only the essential data from one AbstractSolutionData to another.

  Currently copies:
    * q
    * q_vec
    * res
    * res_vec
    * shared_data

  Implementation Notes:
    Avoids double copying when q and q_vec alias
"""
function copyForMultistage!(eqn_dest::AbstractSolutionData, eqn_src::AbstractSolutionData)

  copy!(eqn_dest.q, eqn_src.q)

  # slight optimization- checking if we need to copy q_vec at all,
  #   in case q_vec & q point to the same thing
  if pointer(eqn_src.q) != pointer(eqn_src.q_vec)
    copy!(eqn_dest.q_vec, eqn_src.q_vec)
  end


  # This is needed because we suspect copy! will only work down one level
  num_inner_arrays = length(eqn_src.shared_data)
  for i = 1:length(eqn_src.shared_data)
    copy!(eqn_dest.shared_data[i], eqn_src.shared_data[i])
  end

  copy!(eqn_dest.res, eqn_src.res)
  # same slight optimization as above for q & q_vec
  if pointer(eqn_src.res) != pointer(eqn_src.res_vec)
    copy!(eqn_dest.res_vec, eqn_src.res_vec)
  end

  return nothing
  
end

@doc """

Not in place copy function for AbstractSolutionData. 
    We want this to throw an error because copy! is better.

"""->
function copy(eqn_src::AbstractSolutionData)

  throw(ErrorException("Use copy!(eqn_dest, eqn_src) for AbstractSolution data, not copy"))

  return nothing

end

@doc """
### ODLCommonTools.Boundary

Used to identify boundary faces in a finite-element grid.

**Fields**

* `element` : index of the element to which the boundary face belongs
* `face` : the face index of the boundary (local index to the element)

**Example**

To mark face 2 of element 7 to be a boundary face, use `Boundary(7,2)`

"""->
struct Boundary
  element::UInt32
  face::UInt8
end

"""
  Boundary(0, 0).  Useful as a default value for function arguments
"""
global const NullBoundary = Boundary(0,0)

@doc """
### ODLCommonTools.Interface

Used to identify interfaces between elements in a finite-element grid.

**Fields**

* `elementL` : index of the so-called left element in the pair
* `elementR` : index of the so-called right element in the pair
* `faceL` : the face index of the interface with respect to the left element
* `faceR` : the face index of the interface with respect to the right element
* `orient` : orientation of the 'right' element relative to the 'left'

**Example**

Consider an interface between elements 2 and 5.  Suppose the interface is on
face 1 of element 2 and face 3 of element 5.  Furthermore, suppose element 5 has
orientation 1 relative to element 1 (defintion of orientation TBD).  This can be
indicated as `Interface(2,5,1,3,1)`

"""->
struct Interface
  elementL::UInt32
  elementR::UInt32
  faceL::UInt8
  faceR::UInt8
  orient::UInt8
end

global const NullInterface = Interface(0, 0, 0, 0, 0)

# small interface for Boundary and Interface

"""
  This function returns either the `element` field of a Boundary or the
  `elementL` field of an interface.
"""
function getElementL(bndry::Boundary)
  return bndry.element
end

function getElementL(iface::Interface)
  return iface.elementL
end

"""
  This function returns either the `face` field of a Boundary or the
  `faceL` field of an Interface
"""
function getFaceL(bndry::Boundary)
  return bndry.face
end

function getFaceL(iface::Interface)
  return iface.faceL
end

"""
  Abstract supertype of all boundary condition functors
"""
abstract type BCType end


"""
  Equality test for two BCTypes of the same type.  Tests that all fields
  are equal
"""
function (==)(x::T, y::T) where {T <: BCType}

  fnames = fieldnames(T)
  val = true
  for fname in fnames
    val = val && (getfield(x, fname) == getfield(y, fname))
  end

  return val
end



"""
  Abstract supertype of all boundary condition functors that compute the
  reverse mode with respect to the metrics
"""
abstract type BCType_revm end


"""
  Abstract supertype of all boundary condition functors that compute the
  reverse mode with respect to the solution
"""
abstract type BCType_revq end


"""
  Abstract supertype of all source term functors
"""
abstract type SRCType end

"""
  Abstract supertype of all 2 point numerical flux functions.
"""
abstract type FluxType end # functor DG flux abstract type

"""
  Like [`Fluxtype`](@ref) but for functors that compute the reverse mode
  with respect to the metrics.
"""
abstract type FluxType_revm end # functor type for reverse mode of DG interface fluxes w.r.t mesh metrics


"""
  Like [`FluxType`](@ref), but for flux functions that compute the
  jacobian of the flux with respect to the left and right solutions.
"""
abstract type FluxType_diff end

"""
  Like [`FluxType`](@ref) but for functors that compute the reverse mode
  with respect to the solution q
"""
abstract type FluxType_revq end

abstract type FunctionalType end # functor for functional abstract type

"""
  Show method for Boundary objects
"""
function show(io::IO, object::Boundary)
  print(io, "Boundary element, face = ", object.element, ", ", object.face)
end

"""
  Show method for Interface objects
"""
function show(io::IO, obj::Interface)
  print(io, "Interface elementL, elementR, faceL, faceR, orient = ",
        obj.elementL, ", ",obj.elementR, ", ", obj.faceL, ", ", obj.faceR,
        ", ",obj.orient)
end

"""
  Compare boundaries first by element, then by face
"""
function isless(a::Boundary, b::Boundary)
  if a.element < b.element
    return true
  elseif a.element > b.element
    return false
  elseif a.face > b.face     # elements are same
    return true
  else
    return false
  end

end

"""
  Compare Interfaces, first by elementL, then by elementR
"""
function isless(a::Interface, b::Interface)
  if a.elementL < b.elementL
    return true
  elseif a.elementL > b.elementL
    return false
  end

  if a.elementR < b.elementR
    return true
  else
    return false
  end
end

@doc """
### ODLCommonTools.calcDiffElementArea

  This function calculates the equivalent of the differential area of an element
  at a point in the x and y directions?

  Inputs:
    nrm: a normal vector in the parametric coordinate system of the element 
         at the node
    dxidx: the 2x2 matrix dxi/dx at the node

  Inputs/Outputs:
    workvec: [dx, dy]

"""->
function calcDiffElementArea(nrm::AbstractArray{T,1}, 
                            dxidx::AbstractArray{T2,2},
                            workvec::AbstractArray{T3,1}) where {T, T2, T3}
  fill!(workvec, zero(T3))
  for di1 = 1:size(nrm,1)
    for di2 = 1:size(nrm,1)
      workvec[di2] += nrm[di1]*dxidx[di1,di2]
    end
  end
  return norm(workvec)
end

# it would be better if this used @boundscheck
@doc """
### Utils.safe_views

  This bool value controls whether the function named sview refers to 
  view or unsafe_view from the ArrayViews package
"""->
global const safe_views = false
if safe_views
  global const sview = ArrayViews.aview
else
  global const sview = ArrayViews.unsafe_aview
end

"""
  This macro eliminates blocks of code if safe_views is false
"""
macro ifsafeview(ex)

  if safe_views
    return quote
      $(esc(ex))
    end
  else
    return nothing
  end
end  # end macro

@doc """
### ODLCommonTools.functorThatErrors

  This functor will error.

  It is intended to be used to initialize functor fields in eqn to some value 
    so that eqn_deepcopy can operate upon them.

  This should never be called; if it is called, then the field of eqn has not been
    properly defined for its desired usage.

  Inputs: none
  Outputs: none

"""->
mutable struct functorThatErrors <: FluxType
end

function (obj::functorThatErrors)(params::AbstractParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F::AbstractVector{Tres}) where {Tsol, Tres, Tmsh}

  error("Default functor has been called. You have not properly initialized something.")
  return nothing
end

@doc """
### ODLCommonTools.functorThatErrors_revm

  This functor will error.

  It is intended to be used to initialize functor fields in eqn to some value 
    so that eqn_deepcopy can operate upon them.

  This should never be called; if it is called, then the field of eqn has not been
    properly defined for its desired usage.

  Inputs: none
  Outputs: none

"""->
mutable struct functorThatErrors_revm <: FluxType_revm
end

function (obj::functorThatErrors_revm)( params::AbstractParamType,
              uL::AbstractArray{Tsol,1},
              uR::AbstractArray{Tsol,1},
              aux_vars::AbstractVector{Tres},
              nrm::AbstractVector{Tmsh},
              F::AbstractVector{Tres}) where {Tsol, Tres, Tmsh}

  error("Default functor has been called. You have not properly initialized something.")
  return nothing
end


include("types.jl")
include("abstract_functional.jl")
include("eqn_deepcopy.jl")
include("ro_view.jl")
include("io.jl")
include("misc.jl")
include("getAllTypeParams.jl")
include("lapack.jl")
include("null_array.jl")
include("empty_array.jl")
end     # module
