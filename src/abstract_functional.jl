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
  Evaluates a functional.  This function always returns a value of type
  `Topt`, the static parameter of the [`AbstractFunctional`](@ref) being
  evaluated.  
  """
function evalFunctional
end


"""
  Implementation of [`evalFunctional`](@ref). Each new functional
  should implement this function, not `evalFunctional`.
  The final positional argument of this function *must* be
  the `AbstractFunctional` to be evaluated.
"""
function _evalFunctional
end


"""
  Performs reverse-mode differentiation of a functional with respect to the
  metrics
"""
function evalFunctionalDeriv_m
end


"""
  Implementation of [`evalFunctionalDeriv_m`](@ref).  Each new functional
  should implement this function, not `evalFunctionalDeriv_m`
  The final positional argument of this function *must* be
  the `AbstractFunctional` to be evaluated.


"""
function _evalFunctionalDeriv_m
end


"""
  Computes the derivative of a functional with respect to the solution `q`.
  For efficiency, this should use reverse-mode differentiation internally.
"""
function evalFunctionalDeriv_q
end


"""
  Implementation of [`evalFunctionalDeriv_q`](@ref).  Each new functional
  should implement this function, not `evalFunctionalDeriv_q`.
  The final two positional argument of this function *must* be
  the `AbstractFunctional` to be evaluated and the [`Abstract3DArray`](@ref)
  to be overwritten with the solution
"""
function _evalFunctionalDeriv_q
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


#------------------------------------------------------------------------------
# Composite Functional

"""
  Represents a functional that is the sum of other functionals.  All functional
  must have the same parallel data requirement.

  **Constructors**

  ```
  CompositeFunctional(func1, func2)
  ```

  where `func1` and `func2` are `AbstractFunctionals`](@ref), or

  ```
  CompositeFunctional([func1, func2])
  ```
"""
struct CompositeFunctional{Topt} <: AbstractFunctional{Topt}
  funcs::Vector{AbstractFunctional{Topt}}

  function CompositeFunctional{Topt}(funcs::AbstractVector) where {Topt}
    pdata = getParallelData(funcs[1])
    for i=2:length(funcs)
      if getParallelData(funcs[i]) != PARALLEL_DATA_NONE
        @assert getParallelData(funcs[i]) == pdata
      end
    end

    return new(funcs)
  end
end


function CompositeFunctional(funcs::Vararg{AbstractFunctional{Topt}}) where {Topt}

  return CompositeFunctional{Topt}(collect(funcs))
end


function CompositeFunctional(funcs::AbstractVector{AbstractFunctional{Topt}}) where {Topt}

  return CompositeFunctional{Topt}(funcs)
end

function getParallelData(func::CompositeFunctional)
  return getParallelData(func.funcs[1])
end

# In an ideal world, these would use varags to be fully generic, but
# Julia only allows varags on the final argument, so we have to specify the
# exact signature
#function _evalFunctional(args..., func::CompositeFunctional{Topt}; kwargs...)::Topt where {Topt}
function _evalFunctional(mesh::AbstractMesh{Tmsh},
                         sbp::AbstractOperator, eqn::AbstractSolutionData{Tsol}, opts,
                         func::CompositeFunctional{Topt})::Topt where {Tmsh, Tsol, Topt}

  J = 0
  for f in func.funcs
    J += _evalFunctional(mesh, sbp, eqn, opts, f)
  end

  return J
end


function _evalFunctionalDeriv_m(mesh::AbstractDGMesh{Tmsh},
                         sbp::AbstractOperator, eqn::AbstractSolutionData{Tsol}, opts,
                         func::CompositeFunctional{Topt},
                         val_bar::Number=1) where {Tmsh, Tsol, Topt}


  for f in func.funcs
    _evalFunctionalDeriv_m(mesh, sbp, eqn, opts, f, val_bar)
  end

  return nothing
end


function _evalFunctionalDeriv_q(mesh::AbstractDGMesh{Tmsh},
                         sbp::AbstractOperator, eqn::AbstractSolutionData{Tsol}, opts,
                         func::CompositeFunctional{Topt},
                         func_deriv_arr::Abstract3DArray) where {Tmsh, Tsol, Topt}

  func_deriv_arr2 = zeros(func_deriv_arr)
  for f in func.funcs
    _evalFunctionalDeriv_q(mesh, sbp, eqn, opts, f, func_deriv_arr2)
    for i=1:length(func_deriv_arr)
      func_deriv_arr[i] += func_deriv_arr2[i]
    end
    fill!(func_deriv_arr2, 0)
  end

  return nothing
end
