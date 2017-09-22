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

  # type for testing isFieldDefined
  type TestType
    a::Int
    b::Array{Float64, 1}
    function TestType()
      return new()
    end
  end


facts("--- Testing misc.jl ---") do

  # test mat-vec
  A = rand(3,3)
  x = rand(3)
  b = A*x
  b2 = zeros(3)
  smallmatvec!(A, x, b2)
  b3 = smallmatvec(A, x)

  @fact b --> roughly(b2, atol=1e-14)
  @fact b --> roughly(b3, atol=1e-14)

  A = rand(4, 3)
  x = rand(3)
  b = A*x
  b2 = smallmatvec(A, x)

  @fact b2 --> roughly(b, atol=1e-14)

  # test reverse mode
  A = rand(4,3)

  x = rand(3)
  x2 = zeros(Complex128, 3)
  copy!(x2, x)

  b = zeros(4)
  b2 = zeros(Complex128, 4)

  jac = zeros(4, 3)
  jac2 = zeros(4, 3)
  h = 1e-20
  pert = Complex128(0, h)
  for i=1:3
    x2[i] += pert
    fill!(b2, 0)
    smallmatvec!(A, x2, b2)

    for j=1:4
      jac[j, i] = imag(b2[j])/h
    end

    x2[i] -= pert
  end

  for j=1:4
    b[j] = 1  # actuall b_bar
    fill!(x, 0)  # actuall x_bar
    smallmatvec_revv!(A, x, b)

    for i=1:3
      jac2[j, i] = x[i]
    end

    b[j] = 0
  end

  @fact norm(jac - jac2) --> roughly(0.0, atol=1e-13)

  # test reinterper-BLAS interaction

  println("testing BLAS fallbacks")
  for d=2:20  # size of arrays
    A = rand(Complex128, d, d, 1)
    Av = unsafe_view(A, :, :, 1)
    
    x = rand(d, 1)
    xv = unsafe_view(x, :, 1)
#    xv = x[:, 1]

    b = zeros(Complex128, d, 1)
    bv = unsafe_view(b, :, 1)
#    bv = b[:, 1]

    # test unsafe_view
    smallmatvec!(Av, xv, bv)

    b2 = A[:, :, 1]*x[:, 1]

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    # test safe view
    Av2 = view(A, :, :, 1)
    xv2 = view(x, :, 1)
    bv2 = view(b, :, 1)

    smallmatvec!(Av2, xv2, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    # check mixed argument types
    b3 = zeros(b2)
    smallmatvec!(Av2, xv2, b3)

    @fact norm(b3 - b2) --> roughly(0.0, atol=1e-14)

    smallmatvec!(A[:, :, 1], xv2, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    smallmatvec!(ROView(Av), xv, bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    smallmatvec!(ROView(Av2), xv, bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    smallmatvec!(ROView(Av2), ROView(xv), bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    smallmatvec!(ROView(Av2), xv, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    # test transposed mat-vec
    b2 = A[:, :, 1].'*x[:, 1]

    smallmatTvec!(Av, xv, bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    smallmatTvec!(Av2, xv2, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

  end


  # test mat-mat 
  for d=2:20  # size of arrays
    A = rand(Complex128, d, d, 1)
    Av = unsafe_view(A, :, :, 1)
    
    x = rand(d, d, 1)
    xv = unsafe_view(x, :, :, 1)
#    xv = x[:, 1]

    b = zeros(Complex128, d, d, 1)
    bv = unsafe_view(b, :, :, 1)
#    bv = b[:, 1]

    b2 = A[:, :, 1]*x[:, :, 1]

    smallmatmat!(Av, xv, bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    # test safe view
    Av2 = view(A, :, :, 1)
    xv2 = view(x, :, :, 1)
    bv2 = view(b, :, :, 1)

    smallmatmat!(Av2, xv2, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    b3 = zeros(b2)

    smallmatmat!(Av2, xv2, b3)

    @fact norm(b3 - b2) --> roughly(0.0, atol=1e-14)

    smallmatmat!(A[:, :, 1], xv2, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    smallmatmat!(ROView(Av2), xv, bv)

    @fact norm(bv - bv2) --> roughly(0.0, atol=1e-14)

    smallmatmat!(ROView(Av2), xv, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    # test mat-matT
    b2 = A[:, :, 1]*x[:, :, 1].'
    smallmatmatT!(Av, xv, bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    smallmatmatT!(Av2, xv2, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    smallmatmatT!(A[:, :, 1], xv2, bv2)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    b3 = zeros(b2)
    smallmatmatT!(Av2, xv2, b3)

    @fact norm(b3 - b2) --> roughly(0.0, atol=1e-14)

    smallmatmatT!(ROView(Av2), xv, bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    smallmatmatT!(ROView(Av2), xv, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

    # test matT-mat
    b2 = A[:, :, 1].'*x[:, :, 1]

    smallmatTmat!(Av, xv, bv)

    @fact norm(bv - b2) --> roughly(0.0, atol=1e-14)

    smallmatTmat!(Av2, xv2, bv2)

    @fact norm(bv2 - b2) --> roughly(0.0, atol=1e-14)

  end


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
  smallmatmatT!(A, x, b)
  b3 = smallmatmatT(A, x)
  @fact b --> roughly(b2, atol=1e-14)
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


  # test smallmatTmat
  A = [1. 2 3; 4 5 6; 7 8 9]
  x = 2*A
  b = smallmatTmat(A, x)
  b2 = A.'*x
  @fact b --> roughly(b2, atol=1e-13)

  A = rand(4, 3)
  x = rand(4, 5)
  b = smallmatTmat(A, x)
  b2 = A.'*x
  @fact b --> roughly(b2, atol=1e-13)

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

  face_verts = ElementTopology3().face_verts
  edge_verts = [1 2 3 1 2 3;
                2 3 1 4 4 4]
  topo = ElementTopology{3}(face_verts, edge_verts, topo2=ElementTopology2())

  @fact topo.face_edges --> [1 1 2 3; 2 5 6 6; 3 4 5 4]
  @fact topo.face_edges_flipped --> [false false false true;
                                     false false false false;
                                     false true  true  true]

  # test isFieldDefined

  obj = TestType()

  ret = isFieldDefined(obj, :a)
  @fact ret --> true  # bitstypes fields are always defined

  obj.a = 2
  ret = isFieldDefined(obj, :a)
  @fact ret --> true

  ret = isFieldDefined(obj, :a, :b)
  @fact ret --> false

  @fact_throws isFieldDefined(obj, :a, :c)
  @fact_throws isFieldDefined(obj)

  obj.b = rand(3)
  ret = isFieldDefined(obj, :a, :b)
  @fact ret --> true

  # test file name modification
  fname = "abc.dat"
  fname2 = get_parallel_fname(fname, 1)
  @fact fname2 --> "abc_1.dat"

  fname = "abc"
  fname2 = get_parallel_fname(fname, 1)
  @fact fname2 --> "abc_1"

  # test append_path, prepend_path
  pth = ""
  new_entry = "/dir1/dir2"
  new_pth = append_path(pth, new_entry)
  @fact new_pth --> new_entry

  pth = "/dir0"
  new_pth = append_path(pth, new_entry)
  @fact new_pth --> "/dir0:/dir1/dir2"

  pth = ""
  new_pth = prepend_path(pth, new_entry)
  @fact new_pth --> new_entry

  pth = "/dir0"
  new_pth = prepend_path(pth, new_entry)
  @fact new_pth --> "/dir1/dir2:/dir0"

  str = joinpath_ascii(pwd(), "hello")
  @fact str --> joinpath(pwd(), "hello")


end
