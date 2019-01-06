# add some basic wrappers for Lapack, because Julia Base was too lazy to
# expose an interface that allows reusing all the arrays
import Base.LinAlg.LAPACK.liblapack, Base.LinAlg.LAPACK.getrf!, Base.LinAlg.LAPACK.getrs!, Base.LinAlg.BlasInt
import Base.LinAlg.BLAS.@blasfunc

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
  ccall( (@blasfunc($(Symbol("dgetrf_"))), liblapack), Void,
#  ccall(( :dgetrf_64_, liblapack), Void,
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
  ccall((@blasfunc($(Symbol("dgetrs_"))), liblapack), Void,
        (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
         Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}),
        &trans, &n, &size(B,2), A, &max(1,stride(A,2)), ipiv, B, &max(1,stride(B,2)), info)
  return info[]
end

"""
  Wrapper for DLASWP
"""
function laswp!(A::ContiguousArrays{Float64, 1}, k1::BlasInt, k2::BlasInt, 
                ipiv::AbstractVector{BlasInt}, incx::Integer=1)

  @assert length(ipiv) == k2*abs(incx)

  n = 1
  lda = 1  # length(A) ?
#  lda = max(1, stride(A, 2))
  ccall(( @blasfunc($(Symbol("dlaswp_"))), liblapack), Void, (Ptr{BlasInt}, Ptr{Float64}, 
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, 
            Ptr{BlasInt}), &n, A, &lda, &k1, &k2, ipiv, &incx)

  return nothing

end

end  # end @eval


"""
  This function applies the permutation matrix specified by `ipiv` from 
  `getrf`.  Use `laswp!` for the inverse permutation matrix.
"""
function applyIpiv!(ipiv::AbstractVector{I}, x::AbstractVector) where {I <: Integer}

  @ifsafeview @assert length(x) == length(b)

  # ipiv specifies the swaps to perform for the *inverse* permutation
  # matrix.  Loop backwards to get piv, the non-inverse permutations.
  for i=length(x):-1:1
    destval = x[ipiv[i]]
    x[ipiv[i]] = x[i]
    x[i] = destval
  end

  return nothing
end

import Base.SparseArrays.UMFPACK: UmfpackLU, umf_ctrl, umf_info,
                                  UmfpackIndexTypes, umf_nm, UMFPACK_A,
                                  UMFPACK_At, UMFPACK_Aat, umferror


for itype in UmfpackIndexTypes
  sol_r = umf_nm("solve", :Float64, itype)

  @eval begin
    function solve_suitesparse(U::UmfpackLU{Float64, $itype}, b::AbstractVector{Float64}, 
                              typ::Integer, x::AbstractVector{Float64})

      if length(b) != U.m
        throw(DimensionMismatch())
      end

      ret = ccall( ($sol_r, :libumfpack), $itype, ($itype, Ptr{$itype}, Ptr{$itype}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{Float64}, Ptr{Float64}), typ, U.colptr, U.rowval, U.nzval, x, b, U.numeric, umf_ctrl, umf_info)

      umferror(ret)

      return nothing
    end

  end  # end eval begin
end # loop over itype
