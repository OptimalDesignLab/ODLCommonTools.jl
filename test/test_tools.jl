# test the tools

facts("--- Testing misc.jl ---") do

  A = rand(3,3)
  x = rand(3)
  b = A*x
  b2 = zeros(3)
  smallmatvec!(A, x, b2)

  @fact b => roughly(b2, atol=1e-14)

  A = rand(4, 4)
  x = rand(4, 2)
  b = A*x
  b2 = rand(4, 2)
  smallmatmat!(A, x, b2)

  @fact b => roughly(b2, atol=1e-14)

  A = rand(4, 3)
  x = rand(3, 2)
  b = A*x
  b2 = rand(4, 2)
  smallmatmat!(A, x, b2)

  @fact b=> roughly(b2, atol=1e-14)


  A = rand(3,3)
  x = rand(3,3)
  b = zeros(3,3)

  b2 = A*(x.') 
  smallmatmatT!(A, x, b)

  @fact b => roughly(b2, atol=1e-14)

  A = rand(7,7)
  x = rand(4, 7)
  b = zeros(7, 4)
  b2 = A*(x.')
  println("b2 = ", b2)
  smallmatmatT!(A, x, b)

  @fact b=> roughly(b2, atol=1e-14)


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



  q = FIFOQueue{Int}()

  push!(q, 1)
  push!(q, 2)
  push!(q, 3)

  @fact length(q) => 3
  val = pop!(q)
  @fact val => 3
  @fact isempty(q) => false
  @fact front(q) => 2

  resize!(q, 10)
  @fact length(q.s) => 10

  empty!(q)
  @fact isempty(q) => true


  mat_dense = [1 3 0 0; 2 4 6 0; 0 5 7 0; 0 0 0 8]
  sparse_bnds = [1 1 2 4; 2 3 3 4]

  mat = SparseMatrixCSC(sparse_bnds, Float64)
  @time mat = SparseMatrixCSC(sparse_bnds, Float64)
  mat2 = sparse(mat_dense)

  @fact mat.colptr => mat2.colptr
  @fact mat.rowval => mat2.rowval

  @fact mat2[1,1] => 1
  @fact mat2[2,1] => 2
  @fact mat2[1,2] => 3
  @fact mat2[2,2] => 4
  @fact mat2[3,2] => 5
  @fact mat2[2,3] => 6
  @fact mat2[3,3] => 7
  @fact mat2[4,4] => 8

  mat2[3,2] = 9
  @fact mat2[3,2] => 9

  mat2[4,4] = 10
  @fact mat2[4,4] => 10

end
