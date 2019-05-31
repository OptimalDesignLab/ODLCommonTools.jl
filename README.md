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


# Version History

 * v0.1: stable version
 * v0.2: new parallel constructs and read-only sview
 * v0.3: tag before upgrading from Julia 0.4 to Julia 0.6
 * v0.4: tag after upgrading to Julia 0.6
 * v0.5: new AbstractFunctional API
 * v0.6: move sparsity pattern code to PDESolver, add some new array types

[![Build Status](https://travis-ci.org/OptimalDesignLab/ODLCommonTools.jl.svg)](https://travis-ci.org/OptimalDesignLab/ODLCommonTools.jl)
[![codecov.io](https://codecov.io/github/OptimalDesignLab/ODLCommonTools.jl/coverage.svg?branch=master)](https://codecov.io/github/OptimalDesignLab/ODLCommonTools.jl?branch=master)
