using ODLCommonTools
#using Base.Test
using FactCheck

# write your own tests here
@fact 1 --> 1

include("test_tools.jl")

b1 = Boundary(1, 2)
b2 = Boundary(2, 3)

@fact b1 < b2 --> true
@fact b2 > b1 --> true

b2 = Boundary(1, 3)

@fact b1 < b2 --> true
@fact b2 > b1 --> true

b2 = Boundary(1,2)

@fact (b1 < b2) --> false
@fact (b2 > b1) --> false
