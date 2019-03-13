
@testset "NullArray" begin
  a1 = NullArray{Float64}(1)
  a2 = NullArray{Float64}(1, 2)
  a3 = NullArray{Float64}(1, 2, 3)
  a4 = NullArray{Float64}(1, 2, 3, 4)
  a5 = NullArray{Float64}(1, 2, 3, 4, 5)
  a6 = NullArray{Float64}(1, 2, 3, 4, 5, 6)
  a7 = NullArray{Float64}(1, 2, 3, 4, 5, 6, 7)

  @test size(a1) == (1,)
  @test size(a2) == (1, 2)
  @test size(a3) == (1, 2, 3)
  @test size(a4) == (1, 2, 3, 4)
  @test size(a5) == (1, 2, 3, 4, 5)
  @test size(a6) == (1, 2, 3, 4, 5, 6)
  @test size(a7) == (1, 2, 3, 4, 5, 6, 7)

  arrays = [a1, a2, a3, a4, a5, a6, a7]
  len = 1
  for i=1:length(arrays)
    len *= i
    @test length(arrays[i]) == len
    @test_throws ErrorException arrays[i][1]
    @test_throws ErrorException arrays[i][1] = 1
  end

end

@testset "EmptyArray" begin
  a1 = EmptyArray{Float64}(1)
  a2 = EmptyArray{Float64}(1, 2)
  a3 = EmptyArray{Float64}(1, 2, 3)
  a4 = EmptyArray{Float64}(1, 2, 3, 4)
  a5 = EmptyArray{Float64}(1, 2, 3, 4, 5)
  a6 = EmptyArray{Float64}(1, 2, 3, 4, 5, 6)
  a7 = EmptyArray{Float64}(1, 2, 3, 4, 5, 6, 7)

  @test size(a1) == (1,)
  @test size(a2) == (1, 2)
  @test size(a3) == (1, 2, 3)
  @test size(a4) == (1, 2, 3, 4)
  @test size(a5) == (1, 2, 3, 4, 5)
  @test size(a6) == (1, 2, 3, 4, 5, 6)
  @test size(a7) == (1, 2, 3, 4, 5, 6, 7)

  arrays = [a1, a2, a3, a4, a5, a6, a7]
  len = 1
  v = 0
  for i=1:length(arrays)
    a = arrays[i]
    len *= i
    @test length(a) == len

    # make sure these don't error
    v += a[i]
    a[i] = 1
    a[1] += 1
  end

end


