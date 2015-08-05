module PDESolverCommon

import Base.show

export AbstractSolutionData
export AbstractMesh
export Boundary
export Interface
export BCType

abstract AbstractSolutionData{T3}
abstract AbstractMesh{T1}

@doc """
### PDESolverCommon.Boundary

Used to identify boundary faces in a finite-element grid.

**Fields**

* `element` : index of the element to which the boundary face belongs
* `face` : the face index of the boundary (local index to the element)

**Example**

To mark face 2 of element 7 to be a boundary face, use `Boundary(7,2)`

"""->
immutable Boundary
    element::UInt64
    face::UInt8
end

@doc """
### PDESolverCommon.Interface

Used to identify interfaces between elements in a finite-element grid.

**Fields**

* `elementL` : index of the so-called left element in the pair
* `elementR` : index of the so-called right element in the pair
* `faceL` : the face index of the interface with respect to the left element
* `faceR` : the face index of the interface with respect to the right element
* `nodemap` : local mapping that matches nodes from right element to left

**Example**

Consider an interface between elements 2 and 5.  Suppose the interface is on
face 1 of element 2 and face 3 of element 5.  This can be indicated as
`Interface(2,5,1,3)`

"""->
immutable Interface
  elementL::UInt64
  elementR::UInt64
  faceL::UInt8
  faceR::UInt8
  nodemap::Array{UInt8, 1}
end

abstract BCType  # functor boundary condition abstract type

function show(io::IO, object::Boundary)
  print(io, "Boundary element, face = ", object.element, ", ", object.face)
end

function show(io::IO, obj::Interface)
  print(io, "Interface elementL, elementR, faceL, faceR = ", obj.elementL, ", ",
        obj.elementR, ", ", obj.faceL, ", ", obj.faceR)
end

end # module
