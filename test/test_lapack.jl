import Base.LinAlg.BlasInt
using Base.LinAlg.BLAS


@testset "----- Testing Lapack -----" begin

  A = [1.0 2 3; 4 5 6; 8 8 9]
  A2 = copy(A)
  ipiv = zeros(BlasInt, 3)
  info = getrf!(A, ipiv)
  @test ( info )== 0
  b = [1.0, 2, 3]
  b2 = copy(b)
  info = getrs2!('N', A, ipiv, b)
  @test ( info )== 0
  
  x = copy(b)
  x2 = A2\b2
  @test isapprox( norm(x - x2), 0.0) atol=1e-13

  # test getrs2!
  # test A*x against P*L*U*x

  N = 30
  for i=1:10
    A = rand(N, N)
    x = rand(N)
    b = A*x
    b2 = zeros(b)
    ipiv = zeros(BlasInt, N)
    getrf!(A,ipiv)
    copy!(b2, x)
    trmv!('U', 'N', 'N', A, b2)
    trmv!('L', 'N', 'U', A, b2)
    applyIpiv!(ipiv, b2)

    @test isapprox( norm(b2 - b), 0.0) atol=1e-12
  end

  # test suitesparse solve
  A = [1.0 2 3; 4 5 6; 8 8 9]
  As = sparse(A)
  b = Float64[1, 2, 3]

  x = A\b
  lu = factorize(As)
  x2 = zeros(x)
  solve_suitesparse(lu, b, UMFPACK_A, x2)

  @test isapprox( norm(x2 - x), 0.0) atol=1e-13

end
