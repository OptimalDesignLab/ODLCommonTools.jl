
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
  This function copies the data from one array of arrays to another array
  of arrays with the same structure.  It recurses down to the deepest array
  (whose element type is not a subtype of AbstractArray) and calls copy! on it.

  **Inputs**

   * src: an AbstractArray (possibly an array of arrays)

  **Inputs/Outputs**

   * dest: an AbstractArray (must be similar to src)

  Note that this function will not work correctly if any array (including
  the nested ones) has an element type of Any.

  Aliasing: aliasing is allowed, however cycles in the recursion tree are not
            allowed.  This function does not attempt to handle this case, and
            it will likely result in a StackOverflow error.
"""
function copy_array_recursive!{Td <: AbstractArray, Ts <: AbstractArray}(dest::AbstractArray{Td}, src::AbstractArray{Ts})

  @assert length(dest) == length(src)

  for i=1:length(src)
    # if the element type of src[i] is an AbstractArray, call this method
    # again, otherwise call the other method
    copy_array_recursive!(dest[i], src[i])
  end

  return nothing
end

function copy_array_recursive!{Td, Ts}(dest::AbstractArray{Td},
                                        src::AbstractArray{Ts})

  copy!(dest, src)

  return nothing
end

"""
  ODLCommonTools.eqn_deepcopy_fields

  This function is physics-agnostic and copies over all fields.
  It assumes the fields of the destination object are already initialized
  and (in the case of arrays) the right size.

    **Inputs**:

     * eqn

    **Outputs:**
      
     * eqn_copy

"""
function copy!(eqn_copy::AbstractSolutionData, eqn::AbstractSolutionData)
  # 2: copy over fields

  for fdnm in fieldnames(eqn)    # loop over first level fieldnames in eqn
    println("fieldname = ", fdnm)

    if fdnm == :file_dict
      setfield!(eqn_copy, fdnm, getfield(eqn, fdnm))
    end

    fdnm_type = typeof(getfield(eqn, fdnm))    # get the super type of the current field

    # ------- handle params
    # if this first level fieldname is of type ParamType; ex: eqn.params
    if issubtype(fdnm_type, AbstractParamType)
      params_copy = getfield(eqn_copy, fdnm)
      params = getfield(eqn, fdnm)

      # loop over eqn.params, eqn.params_conservative, or eqn.params_entropy
      # println(" is a subtype of AbstractParamType, fdnm: ", fdnm)

      for fdnm_lvl2 in fieldnames(params)      # loop over 2nd level fieldnames

        # get the super type of the current 2nd level field
        fdnm_lvl2_type = typeof(getfield(params, fdnm_lvl2))

        # if the 2nd level fieldname is of type Array; ex: eqn.params.q_vals
        if issubtype(fdnm_lvl2_type, AbstractArray)

          # this does not work: setfield!(getfield(eqn_copy, a), b , getfield(getfield(eqn, a),b))
          #   must use copy, or else changing eqn's value changes eqn_copy
          copy!(getfield(params_copy, fdnm_lvl2), getfield(params, fdnm_lvl2))
#          setfield!(params_copy, fdnm_lvl2, copy(getfield(params, fdnm_lvl2)))

          # Note: this is assuming that there are no Arrays of Arrays inside 
          #       an eqn.params (or params_entropy, etc)

        else            # if the 2nd level fieldname is not of type Array; ex: eqn.params.gamma

          # because copy is not defined for all non-array types, such as functions
          setfield!(params_copy, fdnm_lvl2, getfield(params, fdnm_lvl2))

        end

      end

    # -------- handle arrays
    # if this first level fieldname is of type Array; ex: eqn.q or eqn.q_face_send
    elseif issubtype(fdnm_type, AbstractArray)

      arr_copy = getfield(eqn_copy, fdnm)
      arr = getfield(eqn, fdnm)

      copy_array_recursive!(arr_copy, arr)
#=
      # -------- handle array of arrays
      # if this is an Array of Arrays; ex: eqn.q_face_send
      if issubtype(eltype(arr), AbstractArray)

        # we only support 2 levels of nested arrays
        @assert !(eltype(arr) <: AbstractArray)

        # first copy the outer array
        # copy is required here, as the innermost object is an array
#        setfield!(eqn_copy, fdnm, copy(getfield(eqn, fdnm)))

        for i=1:length(arr)
          copy!(arr_copy[i], arr[i])
        end
#=
        # then loop over array and copy all elements, which are each an array
        for i = 1:length(getfield(eqn, fdnm))

          setindex!(getfield(eqn_copy, fdnm), getindex(getfield(eqn, fdnm), i), i)

        end   # end of loop over elements (each an array) of the 1st level field, which is of type array
=#
      else      # if this a simple Array; ex: eqn.q

        copy!(arr_copy, arr)
#=
        # handle special case of q_vec and res_vec needing to be reshaped in DG case
        if mesh.isDG && (fdnm == :q_vec)

          eqn_copy.q_vec = reshape(eqn_copy.q, mesh.numDof)

        elseif mesh.isDG && (fdnm == :res_vec)

          eqn_copy.res_vec = reshape(eqn_copy.res, mesh.numDof)

        else
          # copy is required here, as the innermost object is an array
          setfield!(eqn_copy, fdnm, copy(getfield(eqn, fdnm)))
        end
=#
      end
=#
    else        # handle non-arrays

      # copy is not defined for many of these non-array types: use assignment
      setfield!(eqn_copy, fdnm, getfield(eqn, fdnm))

    end

  end     # end of loop over first level fieldnames


  return eqn_copy

end
