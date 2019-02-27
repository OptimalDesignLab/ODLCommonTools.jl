using ODLCommonTools
#using Base.Test
#using FactCheck
using Base.Test
using ArrayViews
import ArrayViews.view

function not_isapprox(args...; kwargs...)
  return !isapprox(args...; kwargs...)
end

@testset "----- Testing Interface and Boundary -----" begin
  # write your own tests here
  @test ( 1 )== 1

  b1 = Boundary(1, 2)
  b2 = Boundary(2, 3)

  @test ( b1 < b2 )== true

  @test ( b2 > b1 )== true

  @test ( getElementL(b1) )== 1
  @test ( getElementL(b2) )== 2
  @test ( getFaceL(b1) )== 2

  b2 = Boundary(1, 3)

  @test ( b1 < b2 )== false
  @test ( b2 > b1 )== false

  b2 = Boundary(1,2)

  @test ( (b1 < b2) )== false
  @test ( (b2 > b1) )== false

  i1 = Interface(1, 2, 3, 4, 5)
  i2 = Interface(2, 3, 4, 5, 6)

  @test ( (i1 < i2) )== true
  @test ( (i2 > i1) )== true

  @test ( getElementL(i1) )== 1
  @test ( getFaceL(i1) )== 3

end

include("test_tools.jl")
include("test_copy.jl")
include("test_roview.jl")
include("test_io.jl")
include("test_lapack.jl")
include("test_null_array.jl")

#FactCheck.exitstatus()
