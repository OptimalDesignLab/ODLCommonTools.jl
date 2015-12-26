# Name: ODLCommonTools
# Description: Common Tools for ODL
#   Performs the function of forward declaring abstract types

module ODLCommonTools
include("misc.jl")
include("sparse.jl")
import Base.show

export AbstractSolutionData, AbstractParamType, Abstract3DArray
export AbstractMesh
export Boundary
export Interface
export BCType
export calcNorm
abstract AbstractSolutionData{Tsol, Tres} # Abstract type defnition
abstract AbstractMesh{Tmsh}
abstract AbstractParamType
typealias Abstract3DArray{T} AbstractArray{T, 3}
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
* `nodemap` : local mapping that matches nodes from right element to left

**Example**

Consider an interface between elements 2 and 5.  Suppose the interface is on
face 1 of element 2 and face 3 of element 5.  This can be indicated as
`Interface(2,5,1,3)`

"""->
immutable Interface
  elementL::UInt32
  elementR::UInt32
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

@doc """
###ODLCommonTools.calcNorm

  This function calculates the norm of a vector (of length numDof) using the
    SBP norm.

    Inputs:
      eqn:  an AbstractSolutionData
      res_vec:  vector to calculate the norm of

    Keyword arguments:
      strongres: if res_vec is the residual of the weak form, then
                 strongres=true computes (efficiently) the norm of the strong
                 form residual.  Default false

    Returns:
      val:  norm of solution using SBP norm (Float64)

    There are no restrctions on the datatype of res_vec (ie. it can be complex)

    Aliasing restrictions: none

"""->
function calcNorm{T}(eqn::AbstractSolutionData, res_vec::AbstractArray{T}; strongres=false)
# calculates the norm of a vector using the mass matrix

  val = zero(real(res_vec[1]))

  if !strongres
    for i=1:length(res_vec)
      val += real(res_vec[i])*eqn.M[i]*real(res_vec[i])   # res^T M res
    end
  else  # strongres
    for i=1:length(res_vec)
      val += real(res_vec[i])*eqn.Minv[i]*real(res_vec[i])   # res^T M res
    end
  end


  val = sqrt(val)
  return val
end     # end of calcNorm function

end     # module
