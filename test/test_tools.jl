  mutable struct TestMesh <: AbstractDGMesh{Float64}
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
  mutable struct TestType
    a::Int
    b::Array{Float64, 1}
    function TestType()
      return new()
    end
  end


function normLinf(A::AbstractArray)

  return maximum(abs.(A))
end


@testset "--- Testing misc.jl ---" begin

  # test mat-vec
  A = rand(3,3)
  x = rand(3)
  b = A*x
  b2 = zeros(3)
  smallmatvec!(A, x, b2)
  b3 = smallmatvec(A, x)

  @test isapprox( b, b2) atol=1e-14
  @test isapprox( b, b3) atol=1e-14

  A = rand(4, 3)
  x = rand(3)
  b = A*x
  b2 = smallmatvec(A, x)

  @test isapprox( b2, b) atol=1e-14

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

  @test isapprox( normLinf(jac - jac2), 0.0) atol=1e-13

  # test reinterpert-BLAS interaction

  println("testing BLAS fallbacks")
  for d=2:20  # size of arrays
    A = rand(Complex128, d, d, 1)
    Av = unsafe_aview(A, :, :, 1)
    
    x = rand(d, 1)
    xv = unsafe_aview(x, :, 1)
#    xv = x[:, 1]

    b = zeros(Complex128, d, 1)
    bv = unsafe_aview(b, :, 1)
#    bv = b[:, 1]

    # test unsafe_aview
    smallmatvec!(Av, xv, bv)

    b2 = A[:, :, 1]*x[:, 1]

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    # test safe view
    Av2 = aview(A, :, :, 1)
    xv2 = aview(x, :, 1)
    bv2 = aview(b, :, 1)

    smallmatvec!(Av2, xv2, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

    # check mixed argument types
    b3 = zeros(b2)
    smallmatvec!(Av2, xv2, b3)

    @test isapprox( normLinf(b3 - b2), 0.0) atol=1e-14

    smallmatvec!(A[:, :, 1], xv2, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

    smallmatvec!(ROView(Av), xv, bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    smallmatvec!(ROView(Av2), xv, bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    smallmatvec!(ROView(Av2), ROView(xv), bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    smallmatvec!(ROView(Av2), xv, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

    # test transposed mat-vec
    b2 = A[:, :, 1].'*x[:, 1]

    smallmatTvec!(Av, xv, bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    smallmatTvec!(Av2, xv2, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

  end


  # test mat-mat 
  for d=2:20  # size of arrays
    A = rand(Complex128, d, d, 1)
    Av = unsafe_aview(A, :, :, 1)
    
    x = rand(d, d, 1)
    xv = unsafe_aview(x, :, :, 1)
#    xv = x[:, 1]

    b = zeros(Complex128, d, d, 1)
    bv = unsafe_aview(b, :, :, 1)
#    bv = b[:, 1]

    b2 = A[:, :, 1]*x[:, :, 1]

    smallmatmat!(Av, xv, bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    # test safe aview
    Av2 = aview(A, :, :, 1)
    xv2 = aview(x, :, :, 1)
    bv2 = aview(b, :, :, 1)

    smallmatmat!(Av2, xv2, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

    b3 = zeros(b2)

    smallmatmat!(Av2, xv2, b3)

    @test isapprox( normLinf(b3 - b2), 0.0) atol=1e-14

    smallmatmat!(A[:, :, 1], xv2, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

    smallmatmat!(ROView(Av2), xv, bv)

    @test isapprox( normLinf(bv - bv2), 0.0) atol=1e-14

    smallmatmat!(ROView(Av2), xv, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14



    # test mat-matT
    b2 = A[:, :, 1]*x[:, :, 1].'
    smallmatmatT!(Av, xv, bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    smallmatmatT!(Av2, xv2, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

    smallmatmatT!(A[:, :, 1], xv2, bv2)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    b3 = zeros(b2)
    smallmatmatT!(Av2, xv2, b3)

    @test isapprox( normLinf(b3 - b2), 0.0) atol=1e-14

    smallmatmatT!(ROView(Av2), xv, bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    smallmatmatT!(ROView(Av2), xv, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

    # test matT-mat
    b2 = A[:, :, 1].'*x[:, :, 1]

    smallmatTmat!(Av, xv, bv)

    @test isapprox( normLinf(bv - b2), 0.0) atol=1e-14

    smallmatTmat!(Av2, xv2, bv2)

    @test isapprox( normLinf(bv2 - b2), 0.0) atol=1e-14

  end

  # A is d1 x d2, x = d2 x d3, b = d1 x d3
  for d1=2:20
    for d2=2:20
      for d3=2:20

        # mat-mat
        A = rand(d1, d2)
        x = rand(d2, d3)
        b = A*x
        b2 = rand(d1, d3)
        smallmatmat!(A, x, b2)
        b3 = smallmatmat(A, x)
        b4 = rand(d1, d3)
        b5 = 2*A*x + 3*b4
        smallmatmat_kernel!(A, x, b4, 2, 3)

        @test normLinf(b - b2) < 1e-13
        @test normLinf(b - b3) < 1e-13
        @test normLinf(b4 - b5) < 1e-13

        # mat-matT
        A = rand(d1, d2)
        x = rand(d3, d2)
        b = A*(x.')
        b2 = rand(d1, d3)
        smallmatmatT!(A, x, b2)
        b3 = smallmatmatT(A, x)
        b4 = rand(d1, d3)
        b5 = 2*A*(x.') + 3*b4
        smallmatmatT_kernel!(A, x, b4, 2, 3)

        @test normLinf(b - b2) < 1e-13
        @test normLinf(b - b3) < 1e-13
        @test normLinf(b4 - b5) < 1e-13

        # matT-mat
        A = rand(d2, d1)
        x = rand(d2, d3)
        b = A.'*x
        b2 = rand(d1, d3)
        smallmatTmat!(A, x, b2)
        b3 = smallmatTmat(A, x)
        b4 = rand(d1, d3)
        b5 = 2*A.'*x + 3*b4
        smallmatTmat_kernel!(A, x, b4, 2, 3)

        @test normLinf(b - b2) < 1e-13
        @test normLinf(b - b3) < 1e-13
        @test normLinf(b4 - b5) < 1e-13
      end
    end
  end

  # test mat-vec
  # A = d1 x d2, x = d2, b = d1
  for d1=1:20
    for d2=1:20
      # mat-vec
      A = rand(d1, d2)
      x = rand(d2)
      b = A*x
      b2 = rand(d1)
      smallmatvec!(A, x, b2)
      b3 = smallmatvec(A, x)
      b4 = rand(d1)
      b5 = 2*A*x + 3*b4
      smallmatvec_kernel!(A, x, b4, 2, 3)

      @test normLinf(b - b2) < 1e-13
      @test normLinf(b - b3) < 1e-13
      @test normLinf(b4 - b5) < 1e-13

      # matT-vec
      A = rand(d2, d1)
      x = rand(d2)
      b = A.'*x
      b2 = rand(d1)
      smallmatTvec!(A, x, b2)
      b3 = smallmatTvec(A, x)
      b4 = rand(d1)
      b5 = 2*A.'*x + 3*b4
      smallmatTvec_kernel!(A, x, b4, 2, 3)

      @test normLinf(b - b2) < 1e-13
      @test normLinf(b - b3) < 1e-13
      @test normLinf(b4 - b5) < 1e-13

    end
  end


  branch_name = getBranchName()
  @test  length(branch_name)  > 0

  time_str = getTimeString()
  @test  length(time_str)  > 0

  A = [0.0 0 1; 0 0 0; 1 0 0]
  numz, arr = checkZeroRows(A, eps())
  @test ( numz )== 1
  @test ( arr )== [false, true, false]

  numz, arr = checkZeroColumns(A, eps())
  @test ( numz )== 1
  @test ( arr )== [false, true, false]

  numz, arr = checkIdenticalColumns(A, 1, eps())
  @test ( numz )== 0
  @test ( arr )== [false, false, false]


  q = FIFOQueue{Int}()

  push!(q, 1)
  push!(q, 2)
  push!(q, 3)

  @test ( length(q) )== 3
  @test ( q.head )== 1
  @test ( q.tail )== 3
  val = pop!(q)
  @test ( val )== 1
  @test ( q.head )== 2
  @test ( isempty(q) )== false
  @test ( front(q) )== 2

  resize!(q, 10)
  @test ( length(q.s) )== 10

  empty!(q)
  @test ( isempty(q) )== true

  q2 = FIFOQueue{Int}(size_hint=100)
  for i=1:100
    push!(q2, i)
  end

  @test ( length(q2) )== 100
  for i=1:50
    @test ( pop!(q2) )== i
  end

  push!(q2, 1)
  @test ( length(q2.s) )== 100
  @test ( length(q2) )== 51
  @test ( pop!(q2) )== 51



  # test calcDiffElementArea
  dxidx = [1. 1; 0 1]
  nrm = [1., 0]
  workvec = rand(2)
  calcDiffElementArea(nrm, dxidx, workvec)
  @test ( workvec )== [1.0, 1.0]

  @test ( rmfile("abc.txt") )== nothing

  s = printbacktrace()
#  println("typeof(s) = ", typeof(s))
#  println("fieldnames(s) = ", fieldnames(s))
  @test ( typeof(s) )== String



  A = [1.0 2.0; 3.0 4.0]
  @test ( isSymmetric(A) )== false
  A[2,1] = A[1,2]
  @test ( isSymmetric(A) )== true
  tol = 1e-12
  A[2,1] += tol/2
  @test ( isSymmetric(A, tol) )== true
  A[2,1] += tol
  @test ( isSymmetric(A, tol) )== false

  make_symmetric!(A)
  @test ( isSymmetric(A, eps()) )== true



  # test topology type
  face_verts = [1 1 1 2; 2 2 3 3; 3 4 4 4]
  ElementTopology{3}(face_verts)  # test the assertions didn't fire

  face_verts[3,2] = 3 # duplicate a face
  @test_throws Exception  ElementTopology{3}(face_verts)

  ElementTopology3()
  ElementTopology2()

  face_verts = ElementTopology3().face_verts
  edge_verts = [1 2 3 1 2 3;
                2 3 1 4 4 4]
  topo = ElementTopology{3}(face_verts, edge_verts, topo2=ElementTopology2())

  @test ( topo.face_edges )== [1 1 2 3; 2 5 6 6; 3 4 5 4]
  @test ( topo.face_edges_flipped )== [false false false true;
                                     false false false false;
                                     false true  true  true]

  # test isFieldDefined

  obj = TestType()

  ret = isFieldDefined(obj, :a)
  @test ( ret )== true  # bitstypes fields are always defined

  obj.a = 2
  ret = isFieldDefined(obj, :a)
  @test ( ret )== true

  ret = isFieldDefined(obj, :a, :b)
  @test ( ret )== false

  @test_throws Exception  isFieldDefined(obj, :a, :c)
  @test_throws Exception  isFieldDefined(obj)

  obj.b = rand(3)
  ret = isFieldDefined(obj, :a, :b)
  @test ( ret )== true

  # test file name modification
  fname = "abc.dat"
  fname2 = get_parallel_fname(fname, 1)
  @test ( fname2 )== "abc_1.dat"

  fname = "abc"
  fname2 = get_parallel_fname(fname, 1)
  @test ( fname2 )== "abc_1"

  # test append_path, prepend_path
  pth = ""
  new_entry = "/dir1/dir2"
  new_pth = append_path(pth, new_entry)
  @test ( new_pth )== new_entry

  pth = "/dir0"
  new_pth = append_path(pth, new_entry)
  @test ( new_pth )== "/dir0:/dir1/dir2"

  pth = ""
  new_pth = prepend_path(pth, new_entry)
  @test ( new_pth )== new_entry

  pth = "/dir0"
  new_pth = prepend_path(pth, new_entry)
  @test ( new_pth )== "/dir1/dir2:/dir0"

  str = joinpath_ascii(pwd(), "hello")
  @test ( str )== joinpath(pwd(), "hello")

  # test splitting fnames
  pth = "test.dat"
  fstem, fext = split_fname(pth)
  @test ( fstem )== "test"
  @test ( fext )== ".dat"

  pth = "test"
  fstem, fext = split_fname(pth)
  @test ( fstem )== "test"
  @test ( fext )== ""






end
