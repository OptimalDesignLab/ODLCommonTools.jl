module PDESolverCommon

import Base.show

# package code goes here
export AbstractEquation
export AbstractMesh
export Boundary
export Interface
export BCType


abstract AbstractEquation{T3}
abstract AbstractMesh{T1}



immutable Boundary
    element::UInt64  # Int32
    face::UInt8
end

immutable Interface
  elementL::UInt64  # Int32
  elementR::UInt64  # int32
  faceL::UInt8
  faceR::UInt8
  nodemap::Array{UInt8, 1}
end


abstract BCType  # functor boundary condition abstract type

function show(io::IO, object::Boundary)
  println(io, "Boundary element, face = ", object.element, " , ", object.face)
end

function show(io::IO, obj::Interface)
  println(io, "Interface elementL, elementR, faceL, faceR = ", obj.elementL, ", ", obj.elementR, ", ", obj.faceL, ", ", obj.faceR)

 end


end # module
