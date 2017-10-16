# add some basic wrappers for Lapack, because Julia Base was too lazy to
# expose an interface that allows reusing all the arrays
import Base.LinAlg.LAPACK.liblapack, Base.LinAlg.LAPACK.getrf!, Base.LinAlg.LAPACK.getrs!, Base.LinAlg.BlasInt
import Base.blasfunc

@eval begin
"""
  A better wrapper for Lapack DGETRF.  The Base function getrs! is adequate
  (although it does error checking every time) for solving the system using
  this factorization.  This function accepts Contiguous ArrayViews

  **Inputs/Outputs**

   * A: matrix to factorize (inplace)
   * ipiv: array of pivot indices

  **Outputs**

   * info: the Lapack info field
"""
function getrf!(A::ContiguousArrays{Float64}, ipiv::AbstractArray{BlasInt})

  # TODO: put this in a debug block
  @ifsafeview begin
    @assert stride(A, 1) == 1
    @assert stride(A, 2) == size(A, 1)
    @assert stride(ipiv, 1) == 1
    @assert length(ipiv) == min( size(A, 1), size(A, 2) )
  end

  info = Ref{BlasInt}()
  m, n = size(A)
  lda  = size(A, 1)
  ccall(( $(blasfunc(:dgetrf_)), liblapack), Void,
        (Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
         Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
        &m, &n, A, &lda, ipiv, info)
  return info[]
end


"""
  A wrapper for DGETRS that does less error checking.  This function accepts
  Contiguous ArrayViews

  **Inputs**

   * trans: char
   * A: the matrix, previously factored by DGETRF
   * ipiv: ipiv from DGETRF

"""
function getrs2!(trans::Char, A::ContiguousArrays{Float64}, ipiv::AbstractVector{BlasInt}, B::ContiguousArrays{Float64})

  @ifsafeview begin
    @assert stride(B, 1) == 1
    @assert stride(B, 2) == size(B, 1)
    @assert size(A, 1) == size(A, 2)
    if size(A, 1) != size(B, 1)
        throw(DimensionMismatch("B has leading dimension $(size(B,1)), but needs $n"))
    end
  end

  m, n = size(A)
  nrhs = size(B, 2)
  info = Ref{BlasInt}()
  ccall(($(blasfunc(:dgetrs_)), liblapack), Void,
        (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
         Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}),
        &trans, &n, &size(B,2), A, &max(1,stride(A,2)), ipiv, B, &max(1,stride(B,2)), info)
  return info[]
end

end  # end @eval
