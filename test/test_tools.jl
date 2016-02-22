# test the tools
  type TestData <: AbstractSolutionData
    M::Array{Float64, 1}
    Minv::Array{Float64, 1}
  end


facts("--- Testing misc.jl ---") do

  A = rand(3,3)
  x = rand(3)
  b = A*x
  b2 = zeros(3)
  smallmatvec!(A, x, b2)
  b3 = smallmatvec(A, x)

  @fact b --> roughly(b2, atol=1e-14)
  @fact b --> roughly(b3, atol=1e-14)

  A = rand(4, 4)
  x = rand(4, 2)
  b = A*x
  b2 = rand(4, 2)
  smallmatmat!(A, x, b2)
  b3 = smallmatmat(A, x)

  @fact b --> roughly(b2, atol=1e-14)
  @fact b --> roughly(b3, atol=1e-14)

  A = rand(4, 3)
  x = rand(3, 2)
  b = A*x
  b2 = rand(4, 2)
  smallmatmat!(A, x, b2)
  b3 = smallmatmat(A, x)
  @fact b--> roughly(b2, atol=1e-14)
  @fact b --> roughly(b3, atol=1e-14)

  A = rand(3,3)
  x = rand(3,3)
  b = zeros(3,3)

  b2 = A*(x.') 
  smallmatmatT!(A, x, b)

  @fact b --> roughly(b2, atol=1e-14)

  A = rand(7,7)
  x = rand(4, 7)
  b = zeros(7, 4)
  b2 = A*(x.')
  println("b2 = ", b2)
  smallmatmatT!(A, x, b)
  b3 = smallmatmatT(A, x)
  @fact b--> roughly(b2, atol=1e-14)
  @fact b --> roughly(b3, atol=1e-14)

  # test smallmatTvec
  A = rand(3,3)
  x = rand(3)
  b = smallmatTvec(A, x)
  b2 = A.'*x
  @fact b --> roughly(b2, atol=1e-14)

  A = rand(4, 3)
  x = rand(4)
  b2 = A.'*x
  b = smallmatTvec(A, x)
  @fact b --> roughly(b2, atol=1e-14)


  A = [0.0 0 1; 0 0 0; 1 0 0]
  numz, arr = checkZeroRows(A, eps())
  @fact numz --> 1
  @fact arr --> [false, true, false]

  numz, arr = checkZeroColumns(A, eps())
  @fact numz --> 1
  @fact arr --> [false, true, false]

  numz, arr = checkIdenticalColumns(A, 1, eps())
  @fact numz --> 0
  @fact arr --> [false, false, false]



  q = FIFOQueue{Int}()

  push!(q, 1)
  push!(q, 2)
  push!(q, 3)

  @fact length(q) --> 3
  val = pop!(q)
  @fact val --> 3
  @fact isempty(q) --> false
  @fact front(q) --> 2

  resize!(q, 10)
  @fact length(q.s) --> 10

  empty!(q)
  @fact isempty(q) --> true


  mat_dense = [1 3 0 0; 2 4 6 0; 0 5 7 0; 0 0 0 8]
  sparse_bnds = [1 1 2 4; 2 3 3 4]

  mat = SparseMatrixCSC(sparse_bnds, Float64)
  @time mat = SparseMatrixCSC(sparse_bnds, Float64)
  mat2 = sparse(mat_dense)

  @fact mat.colptr --> mat2.colptr
  @fact mat.rowval --> mat2.rowval

  @fact mat2[1,1] --> 1
  @fact mat2[2,1] --> 2
  @fact mat2[1,2] --> 3
  @fact mat2[2,2] --> 4
  @fact mat2[3,2] --> 5
  @fact mat2[2,3] --> 6
  @fact mat2[3,3] --> 7
  @fact mat2[4,4] --> 8

  mat2[3,2] = 9
  @fact mat2[3,2] --> 9

  mat2[4,4] = 10
  @fact mat2[4,4] --> 10

  # functions for checking sparsity
  cnt, out_of_bounds = checkSparseColumns(mat_dense, sparse_bnds, 1e-14)
  @fact cnt --> 0

  cnt, out_of_bounds = checkSparseRows(mat_dense, sparse_bnds, 1e-14)
  @fact cnt --> 0

  # modify sparse_bnds to verify the functions detect out of bounds entries
  sparse_bnds2 = copy(sparse_bnds)
  sparse_bnds2[2, 1] = 1
  cnt, out_of_bounds = checkSparseColumns(mat_dense, sparse_bnds2, 1e-14)
  @fact cnt --> 1
  @fact out_of_bounds --> [true, false, false, false]

  cnt, out_of_bounds = checkSparseRows(mat_dense, sparse_bnds2, 1e-14)
  @fact cnt --> 1
  @fact out_of_bounds --> [true, false, false, false]


  cnt = findLarge(mat_dense, 7.5)
  @fact cnt --> 1




  M = rand(10)
  Minv = 1./M
  obj = TestData(M, Minv)
  data = rand(10)

  norm1 = calcNorm(obj, data)
  norm2 = sqrt(sum(data.*M.*data))

  @fact norm1 --> roughly(norm2)

  norm1 = calcNorm(obj, data, strongres=true)
  norm2 = sqrt(sum(data.*Minv.*data))


  data = complex(rand(10), rand(10))
  norm1 = calcNorm(obj, data)
  norm2 = sqrt(sum(real(data).*M.*real(data)))

  @fact norm1 --> roughly(norm2)

  norm1 = calcNorm(obj, data, strongres=true)
  norm2 = sqrt(sum(real(data).*Minv.*real(data)))

  @fact norm1 --> roughly(norm2)

  # test calcDiffElementArea
  dxidx = [1. 1; 0 1]
  nrm = [1., 0]
  workvec = rand(2)
  calcDiffElementArea(nrm, dxidx, workvec)
  @fact workvec --> [1.0, 1.0]

  @fact rmfile("abc.txt") --> nothing

  s = printbacktrace()
#  println("typeof(s) = ", typeof(s))
#  println("fieldnames(s) = ", fieldnames(s))
  @fact typeof(s) --> ASCIIString



  A = [1.0 2.0; 3.0 4.0]
  @fact isSymmetric(A) --> false
  A[2,1] = A[1,2]
  @fact isSymmetric(A) --> true
  tol = 1e-12
  A[2,1] += tol/2
  @fact isSymmetric(A, tol) --> true
  A[2,1] += tol
  @fact isSymmetric(A, tol) --> false

  make_symmetric!(A)
  @fact isSymmetric(A, eps()) --> true
end
