
#=
function eqn_deepcopy(eqn::AbstractSolutionData)

  error("Generic fallback eqn_deepcopy() called erroneously. eqn_deepcopy() should be extended for each physics module's corresponding eqn object.")

  return nothing

end
=#

"""
  ODLCommonTools.eqn_deepcopy

  This function performs a proper deepcopy (unlike julia's builtin deepcopy) 
    on an equation object.
  It preserves reference topology (i.e. q & q_vec pointing to same array in DG schemes).

    Inputs:
      eqn
      mesh
      sbp
      opts

    Outputs:
      eqn_copy

    One reason for doing this is this case:
      a = rand(2,2)
      b = a
      a[3] = 8
      b[3] == 8
      this is because 'a[3] =' is actually setindex!

"""
# function eqn_deepcopy{Tmsh, Tsol, Tres}(mesh::AbstractMesh{Tmsh}, sbp, eqn::AbstractSolutionData{Tsol, Tres}, opts::Dict)
function eqn_deepcopy(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData, opts::Dict)

  # The eqn object has over 100 fields, so it is necessary to write a better approach for 
  #   copying than explicitly copying every named field
  # This is the second write of eqn_deepcopy; the first version explicitly copied each field. 
  #   Was done for SimpleODE and Advection before Euler made it obvious that that was intractable.

  error("eqn_deepcopy fallback reached. A physics-specific method should be defined.")
 
  return nothing

end

"""
  ODLCommonTools.eqn_deepcopy_fields

  This function is physics-agnostic and copies over all fields.

    Inputs:
      mesh
      sbp
      eqn
      eqn_copy
      opts

    Outputs:
      eqn_copy

"""
function eqn_deepcopy_fields(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                             eqn_copy::AbstractSolutionData, opts::Dict)
  # 2: copy over fields

  for fdnm in fieldnames(eqn)    # loop over first level fieldnames in eqn

    fdnm_type = typeof(getfield(eqn, fdnm))    # get the super type of the current field

    # ------- handle params
    # if this first level fieldname is of type ParamType; ex: eqn.params
    if issubtype(fdnm_type, AbstractParamType)

      # loop over eqn.params, eqn.params_conservative, or eqn.params_entropy
      # println(" is a subtype of AbstractParamType, fdnm: ", fdnm)

      for fdnm_lvl2 in fieldnames(getfield(eqn, fdnm))      # loop over 2nd level fieldnames

        # get the super type of the current 2nd level field
        fdnm_lvl2_type = typeof(getfield(getfield(eqn, fdnm), fdnm_lvl2))

        # if the 2nd level fieldname is of type Array; ex: eqn.params.q_vals
        if issubtype(fdnm_lvl2_type, AbstractArray)

          # this does not work: setfield!(getfield(eqn_copy, a), b , getfield(getfield(eqn, a),b))
          #   must use copy, or else changing eqn's value changes eqn_copy
          setfield!(getfield(eqn_copy, fdnm), fdnm_lvl2, copy(getfield(getfield(eqn, fdnm), fdnm_lvl2)))

          # Note: this is assuming that there are no Arrays of Arrays inside 
          #       an eqn.params (or params_entropy, etc)

        else            # if the 2nd level fieldname is not of type Array; ex: eqn.params.gamma

          # Debug statements:
          # println("fdnm: ", fdnm)
          # println("fdnm_lvl2: ", fdnm_lvl2)

          # because copy is not defined for all non-array types, such as functions
          setfield!(getfield(eqn_copy, fdnm), fdnm_lvl2, getfield(getfield(eqn, fdnm), fdnm_lvl2))

        end

      end

    # -------- handle arrays
    # if this first level fieldname is of type Array; ex: eqn.q or eqn.q_face_send
    elseif issubtype(fdnm_type, AbstractArray)

      # -------- handle array of arrays
      # if this is an Array of Arrays; ex: eqn.q_face_send
      if issubtype(eltype(fdnm_type), AbstractArray)

        # first copy the outer array
        # copy is required here, as the innermost object is an array
        setfield!(eqn_copy, fdnm, copy(getfield(eqn, fdnm)))

        # then loop over array and copy all elements, which are each an array
        for i = 1:length(getfield(eqn, fdnm))

          setindex!(getfield(eqn_copy, fdnm), getindex(getfield(eqn, fdnm), i), i)

        end   # end of loop over elements (each an array) of the 1st level field, which is of type array

      else      # if this a simple Array; ex: eqn.q

        # handle special case of q_vec and res_vec needing to be reshaped in DG case
        if mesh.isDG && (fdnm == :q_vec)

          eqn_copy.q_vec = reshape(eqn_copy.q, mesh.numDof)

        elseif mesh.isDG && (fdnm == :res_vec)

          eqn_copy.res_vec = reshape(eqn_copy.res, mesh.numDof)

        else
          # copy is required here, as the innermost object is an array
          setfield!(eqn_copy, fdnm, copy(getfield(eqn, fdnm)))
        end

      end

    else        # handle non-arrays

      # copy is not defined for many of these non-array types: use assignment
      setfield!(eqn_copy, fdnm, getfield(eqn, fdnm))

    end

  end     # end of loop over first level fieldnames


  return eqn_copy

end
