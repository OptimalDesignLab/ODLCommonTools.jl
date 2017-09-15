@doc """
### ODLCommonTools.getAllTypeParams

Default fallback method for getAllTypeParams.
  
  Input: 
    mesh
    eqn
    opts

  Output:
    nothing
"""->
function getAllTypeParams(mesh::AbstractMesh, eqn::AbstractSolutionData, opts)

  error("getAllTypeParams default fallback reached. A physics-specific definition for this is needed.")

  return nothing

end
