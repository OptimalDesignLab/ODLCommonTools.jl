  type TestMesh <: AbstractDGMesh
    pertNeighborEls::Array{Int, 2}
    dofs::Array{Int, 3}
    neighbor_nums::Array{Int, 2}
    numDof::Int
    numNodesPerElement::Int
    numDofPerNode::Int
    numEl::Int
    coloringDistance::Int
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

  branch_name = getBranchName()
  @fact length(branch_name) --> greater_than(0)

  time_str = getTimeString()
  @fact length(time_str) --> greater_than(0)

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
  @fact q.head --> 1
  @fact q.tail --> 3
  val = pop!(q)
  @fact val --> 1
  @fact q.head --> 2
  @fact isempty(q) --> false
  @fact front(q) --> 2

  resize!(q, 10)
  @fact length(q.s) --> 10

  empty!(q)
  @fact isempty(q) --> true

  q2 = FIFOQueue{Int}(size_hint=100)
  for i=1:100
    push!(q2, i)
  end

  @fact length(q2) --> 100
  for i=1:50
    @fact pop!(q2) --> i
  end

  push!(q2, 1)
  @fact length(q2.s) --> 100
  @fact length(q2) --> 51
  @fact pop!(q2) --> 51


  # test SparseMatrixCSC
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



  # check tighter SparseMatrixCSC
  coloringDistance = 2
  pertNeighborEls = [1 0; 2 0]
  neighbor_nums = Array(Int, 1, 1)
  numNodesPerElement = 3
  numDofPerNode = 4
  numEl = 2
  numDof = numDofPerNode*numNodesPerElement*numEl
  dofs = zeros(Int, numDofPerNode, numNodesPerElement, numEl)
  for i=1:length(dofs)
    dofs[i] = i
  end


  mesh = TestMesh(pertNeighborEls, dofs, neighbor_nums, numDof, numNodesPerElement, numDofPerNode, numEl, coloringDistance)

  mat = SparseMatrixCSC(mesh, Int)

  # verify sortedness of rows
  for i=1:mesh.numDof
    start_ptr = mat.colptr[i]
    end_ptr = mat.colptr[i+1] - 1
    vals_i = mat.rowval[start_ptr:end_ptr]
    @fact issorted(vals_i) --> true
  end

  mat_dense = zeros(numDof, numDof)
  mat_dense[1:12, 1:12] = 1
  mat_dense[13:24, 13:24] = 2
  mat_sparse = sparse(mat_dense)
  println("mat.colptr = ", mat.colptr)
  println("mat_sparse.colptr = ", mat_sparse.colptr)
  println("length (mat.rowval) = ", length(mat.rowval))
  println("length(mat_sparse.rowval) = ", length(mat_sparse.rowval))

  println("mat.rowval = ", mat.rowval)
  println("mat_sparse.rowval = ", mat_sparse.rowval)
  @fact mat.colptr --> mat_sparse.colptr
  @fact mat.rowval --> mat_sparse.rowval





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


  # test fastfind
  A = rand(Int, 999)
  sort!(A)
  for i=1:length(A)
    idx = ODLCommonTools.fastfind(A, A[i])
    @fact idx --> i
  end

  A = rand(Int, 1000)
  sort!(A)
  for i=1:length(A)
    idx = ODLCommonTools.fastfind(A, A[i])
    @fact idx --> i
  end


  A = collect(1:10)
  idx = ODLCommonTools.fastfind(A, 11)
  @fact idx --> 0

  A = collect(1:2)
  idx = ODLCommonTools.fastfind(A, 2)
  @fact idx --> 2


  # test topology type
  face_verts = [1 1 1 2; 2 2 3 3; 3 4 4 4]
  ElementTopology{3}(face_verts)  # test the assertions didn't fire

  face_verts[3,2] = 3 # duplicate a face
  @fact_throws ElementTopology{3}(face_verts)

  ElementTopology3()
  ElementTopology2()

end
