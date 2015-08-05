# test the tools

facts("--- Testing misc.jl ---") do

  A = rand(3,3)
  x = rand(3)
  b = A*x
  b2 = zeros(3)
  smallmatvec!(A, x, b2)

  @fact b => roughly(b2, atol=1e-14)

  A = [0.0 0 1; 0 0 0; 1 0 0]
  numz, arr = checkZeroRows(A, eps())
  @fact numz => 1
  @fact arr => [false, true, false]

  numz, arr = checkZeroColumns(A, eps())
  @fact numz => 1
  @fact arr => [false, true, false]

  numz, arr = checkIdenticalColumns(A, 1, eps())
  @fact numz => 0
  @fact arr => [false, false, false]



end
