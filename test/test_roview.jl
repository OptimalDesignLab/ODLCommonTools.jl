# test ro_view
#using ODLCommonTools
#using FactCheck

function sum_cols(A::AbstractMatrix{T}) where T

  val = zero(T)
  for i=1:size(A, 2)
    tmp = ro_sview(A, :, i)
    val += sum(tmp)
  end

  return val
end

function sumit(A::AbstractArray{T}) where T

  val = zero(T)
  for i=1:length(A)
    val += A[i]
  end

  return val
end

@testset "----- Testing ROView -----" begin

  a = rand(100, 1000)

  b = ROView(a)
  @test ( size(b) )== size(a)
  @test ( length(b) )== length(a)

  @test isapprox( norm(a - b), 0.0) 

  # check that allocating a ROView does not allocate 
  if !ODLCommonTools.safe_views
    @time sum_cols(a)
    bytes_alloc = @allocated sum_cols(a)

    @test  bytes_alloc  < 100
  end

  # check that indexing a ROView does not allocate
  if Base.JLOptions().can_inline == true
    sumit(b)
    bytes_alloc = @allocated sumit(b)
    @test  bytes_alloc  < 100
  end

  # test reinterpert
  A = rand(9, 9)
  Av = unsafe_aview(A, 6:9, 1)

  Av2 = reinterpret(Complex128, Av, (2, 1))

  @test ( Av2[1] )== Complex128(Av[1], Av[2])
  @test ( Av2[2] )== Complex128(Av[3], Av[4])

end

