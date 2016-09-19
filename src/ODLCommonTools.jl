__precompile__(true)
# Name: ODLCommonTools
# Description: Common Tools for ODL
#   Performs the function of forward declaring abstract types

module ODLCommonTools
import Base.show
import Base.isless
export AbstractSolutionData, AbstractParamType, Abstract3DArray
export AbstractMesh, AbstractCGMesh, AbstractDGMesh
export Boundary
export Interface
export BCType, SRCType, FluxType, FunctionalType
export calcNorm, calcDiffElementArea
export ElementTopology3, ElementTopology2, ElementTopology

@doc """
### ODLCommonTools.AbtractSolutionData{Tsol, Tres}

  This abstract type is the supertype for all the objects that store the 
  solution data. Every physics module should implement its own subtype.

  Static parameters:
    Tsol: datatype of solution variables
    Tres: datatype of the mesh variables

  See the repo_root/doc/interfaces.md for the description of everything this
  type must implement.

"""->
abstract AbstractSolutionData{Tsol, Tres} # Abstract type defnition
@doc """
### ODLCommonTools.AbstractMesh{Tmsh}

  This abstract type is the supertype for all mesh objects.  Every interface to
  a mesh software should define its own implementation.

  Static parameters:
    Tmsh: datatype of the mesh data (coordinates, mapping to/from parametric
          space, mapping jacobian).

  See the repo_root/doc/interfaces.md for the description of everything this
  type must implement.

"""->
abstract AbstractMesh{Tmsh}

@doc """
### ODLCommonTools.AbstractCGMesh

  The abstrac type is the supertype of all continuous Galerkin meshes
"""->
abstract AbstractCGMesh{Tmsh} <: AbstractMesh{Tmsh}

@doc """
### ODLCommonTools.AbstractDGGMesh

  The abstrac type is the supertype of all discontinuous Galerkin meshes
"""->
abstract AbstractDGMesh{Tmsh} <: AbstractMesh{Tmsh}

@doc """
### ODLCommonTools.AbstractParamType

  This abstract type is the supertype for all Param objects, which hold values 
  needed for the computation in a place that is fast to access.

   The Param type is also useful for dispatching to low level functions which 
   the AbstractSolutionData might not be passed (depending on the organization 
   of the physics module.

"""->
abstract AbstractParamType

@doc """
### ODLCommonTools.Abstract3DArray

  A typealias useful for specify a 3 dimensional AbstractArray without 
  specifying the element datatype.

"""->
typealias Abstract3DArray{T} AbstractArray{T, 3}

@doc """
### ODLCommonTools.AbstractOptimizationData

Abstract datatype for optimization related data structures. All data types
corresponding to optimization problems should be a subtype of this.

"""->
abstract AbstractOptimizationData


@doc """
### ODLCommonTools.Boundary

Used to identify boundary faces in a finite-element grid.

**Fields**

* `element` : index of the element to which the boundary face belongs
* `face` : the face index of the boundary (local index to the element)

**Example**

To mark face 2 of element 7 to be a boundary face, use `Boundary(7,2)`

"""->
immutable Boundary
  element::UInt32
  face::UInt8
end

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
immutable Interface
  elementL::UInt32
  elementR::UInt32
  faceL::UInt8
  faceR::UInt8
  orient::UInt8
end

abstract BCType  # functor boundary condition abstract type

abstract SRCType # functor source term abstract type

abstract FluxType # functor DG flux abstract type

abstract FunctionalType # functor for functional abstract type

function show(io::IO, object::Boundary)
  print(io, "Boundary element, face = ", object.element, ", ", object.face)
end

function show(io::IO, obj::Interface)
  print(io, "Interface elementL, elementR, faceL, faceR, orient = ",
        obj.elementL, ", ",obj.elementR, ", ", obj.faceL, ", ", obj.faceR,
        ", ",obj.orient)
end


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

"""
### ODLCommonTools.ElementTopology3

  This type describes the topology of the reference element.  For now, it only
  contains face information, but eventually will needed more info.
"""

immutable ElementTopology{Tdim}
  face_verts::Array{Int, 2}  # 3 x 4 array holding the vertices of each face of                                a tet

  function ElementTopology(face_verts::Array{Int, 2})

    # do sanity checks

    # check all vertices are within range
    for i=1:length(face_verts)
      @assert face_verts[i] > 0
      @assert face_verts[i] <= 4
    end

    # check faces are distinct
    for i=1:size(face_verts, 2)
      curr_face = sort(face_verts[:, i])
      for j=(i+1):size(face_verts, 2)
        @assert sort(face_verts[:, j]) != curr_face
      end
    end

    return new(face_verts)

  end  # end function
end

"""
  Default constructor that uses Pumi topology
"""
function ElementTopology3()
  face_verts = [1 1 2 1; 2 2 3 3; 3 4 4 4]
  return ElementTopology{3}(face_verts)
end

"""
  Constructs a dummy 2d ElementTopology
"""
function ElementTopology2()
  face_verts = [1 2 3; 2 3 1]
  return ElementTopology{2}(face_verts)
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
function calcDiffElementArea{T, T2, T3}(nrm::AbstractArray{T,1}, 
                                       dxidx::AbstractArray{T2,2},
                                       workvec::AbstractArray{T3,1})
  fill!(workvec, zero(T3))
  for di1 = 1:size(nrm,1)
    for di2 = 1:size(nrm,1)
      workvec[di2] += nrm[di1]*dxidx[di1,di2]
    end
  end
  return norm(workvec)
end


include("misc.jl")
include("sparse.jl")
end     # module
