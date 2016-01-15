# ODLCommonTools
 
This repository contains auxiliary functions useful for all Optimal Design Lab Projects.

## ODLCommonTools.jl
In particular, it contains the abstract type definitions for the PDESolver 
project, as well as auxiliary types and functions, such as `Boundary` 
`Interface`, `calcNorm`, and `calcDiffElementArea`.

## misc.jl
In addition, it contains matrix-vector and matrix-matrix multiplication 
routines optimized for small matrices and bench-marked against the built-in
Julia equivalents (which call the BLAS implementation used by Julia)

## sparse.jl
Several functions for creating and using sparse matrices exist.
Note that these functions override the existing definitions of 
`setindex!` and `getindex`, to versions optimized for variable
bandwidth dense-with-the-band `SparseMatrixCSC`s.


[![Build Status](https://travis-ci.org/OptimalDesignLab/ODLCommonTools.jl.svg)](https://travis-ci.org/OptimalDesignLab/ODLCommonTools.jl)
[![Coverage Status](https://coveralls.io/repos/OptimalDesignLab/ODLCommonTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/OptimalDesignLab/ODLCommonTools.jl?branch=master)

[![Coverage Status](https://coveralls.io/repos/OptimalDesignLab/ODLCommonTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/OptimalDesignLab/ODLCommonTools.jl?branch=master)
[![codecov.io](https://codecov.io/github/OptimalDesignLab/ODLCommonTools.jl/coverage.svg?branch=master)](https://codecov.io/github/OptimalDesignLab/ODLCommonTools.jl?branch=master)
