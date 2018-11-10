# AbstractFunctional

# allow sending of arbitrary face/element-based data
global const PARALLEL_DATA_NONE = 1000
global const PARALLEL_DATA_FACE = 1001
global const PARALLEL_DATA_ELEMENT = 1002


"""
  Given the string describing what data is shared in parallel, returns the
  enum value.  An exception is thrown if an unrecognized value is supplied.

  **Inputs**

   * pdata: a string describing what data is shared in parallel

  **Outputs**

   * an Int enum
"""
function getParallelDataEnum(pdata::String)
  if pdata == "none"
    return PARALLEL_DATA_NONE
  elseif pdata == "face"
    return PARALLEL_DATA_FACE
  elseif pdata == "element"
    return PARALLEL_DATA_ELEMENT
  else
    error("unrecognized parallel data string: $pdata")
  end
end

"""
  Inverse of [`getParallelDataEnum`](@ref), takes the enum and returns the
  string.

  **Inputs**

   * pdata: Int enum

  **Outputs**

   * a String describing the parallel data
"""
function getParallelDataString(pdata::Int)
  if pdata == PARALLEL_DATA_NONE
    return "none"
  elseif pdata == PARALLEL_DATA_FACE
    return "face"
  elseif pdata == PARALLEL_DATA_ELEMENT
    return "element"
  else
    error("unrecognized parallel data enum: $pdata")
  end
end




@doc """
### ODLCommonTools.AbstractFunctional

Abstract datatype for optimization related data structures. All data types
corresponding to optimization problems should be a subtype of this.

  **Static Parameters**

   * Topt

"""->
abstract type AbstractFunctional{Topt} end


"""
  Abstract type for integral objective functions.  Any integral objective
  function should be a subtype of this, which is a subtype of
  [`AbstractFunctional`](@ref).

  **Static Parameters**

   * Topt
"""
abstract type AbstractIntegralFunctional{Topt} <: AbstractFunctional{Topt} end

"""
  Abstract type for functionals that are boundary integrals

  **Static Parameters**

   * Topt
"""
abstract type AbstractBoundaryFunctional{Topt} <: AbstractIntegralFunctional{Topt} end



#------------------------------------------------------------------------------
# API for functionals

# evalFunctional, evalFunctionalDeriv_q and evalFunctionalDeriv_m are defined
# in more detail in PDESolver

"""
  Evaluates a functional
"""
function evalFunctional
end


"""
  Performs reverse-mode differentiation of a functional with respect to the
  metrics
"""
function evalFunctionalDeriv_m
end


"""
  Computes the derivative of a functional with respect to the solution `q`.
  For efficiency, this should use reverse-mode differentiation internally.
"""
function evalFunctionalDeriv_q
end


"""
  A setup function called immediately before the functional or its derivatives
  are calculated.  This provides the opportunity to cache any
  expensive-to-compute values inside the functional object that might be
  needed when evaluating the integrand.

  **Inputs**

   * mesh
   * sbp
   * eqn
   * opts
"""
function setupFunctional(mesh::AbstractMesh, sbp,
                         eqn::AbstractSolutionData, opts::Dict,
                         functional::AbstractFunctional)

  # in general nothing to do here
  return nothing
end


"""
  Function that describes what kind of solution parallel data is required to
  evaluate this functional.  The function should return one of:

   * `PARALLEL_DATA_NONE`
   * `PARALLEL_DATA_FACE`
   * `PARALLEL_DATA_ELEMENT`

  Note that it is not necessary to implement this method for
  [`AbstractBoundaryFunctional`](@ref).

  **Inputs**

   * obj: the [`AbstractFunctional`](@ref)

  **Outputs**

   * one of the above strings
"""
function getParallelData(obj::AbstractFunctional)

  error("generic fallback for getParallelData(::AbstractFunctional) reached.  Did you forget to extend getParallelData with a new method for your functional?")

end


function getParallelData(obj::AbstractBoundaryFunctional)

  return PARALLEL_DATA_NONE
end
