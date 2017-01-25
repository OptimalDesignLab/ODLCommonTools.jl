using ODLCommonTools
#using Base.Test
using FactCheck

# write your own tests here
@fact 1 --> 1

include("test_tools.jl")
include("test_copy.jl")

b1 = Boundary(1, 2)
b2 = Boundary(2, 3)

@fact b1 < b2 --> true
@fact b2 > b1 --> true

@fact getElementL(b1) --> 1
@fact getElementL(b2) --> 2
@fact getFaceL(b1) --> 2

b2 = Boundary(1, 3)

@fact b1 < b2 --> true
@fact b2 > b1 --> true

b2 = Boundary(1,2)

@fact (b1 < b2) --> false
@fact (b2 > b1) --> false

i1 = Interface(1, 2, 3, 4, 5)
i2 = Interface(2, 3, 4, 5, 6)

@fact (i1 < i2) --> true
@fact (i2 > i1) --> true

@fact getElementL(i1) --> 1
@fact getFaceL(i1) --> 3
