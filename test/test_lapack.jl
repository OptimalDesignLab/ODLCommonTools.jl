import Base.LinAlg.BlasInt


facts("----- Testing Lapack -----") do

  A = [1.0 2 3; 4 5 6; 8 8 9]
  A2 = copy(A)
  ipiv = zeros(BlasInt, 3)
  info = getrf!(A, ipiv)
  @fact info --> 0
  b = [1.0, 2, 3]
  b2 = copy(b)
  info = getrs2!('N', A, ipiv, b)
  @fact info --> 0
  
  x = copy(b)
  x2 = A2\b2
  @fact norm(x - x2) --> roughly(0.0, atol=1e-13)

end
