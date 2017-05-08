# eqn_copy.jl
#
# Contains the generic fallback methods for eqn object manipulation:
#   copy()
#   getStaticParams()
#
# The methods defined here are intended to be extended by physics modules,
#   and therefore should only error out here.
#

function copy(eqn::AbstractSolutionData)

  error("Generic fallback copy() called erroneously. copy() should be extended for each physics module's corresponding eqn object.")

  return nothing

end

function getStaticParams(eqn::AbstractSolutionData)

  error("Generic fallback getStaticParams() called erroneously. getStaticParams() should be extended for each physics module's corresponding eqn object.")

  return nothing

end


