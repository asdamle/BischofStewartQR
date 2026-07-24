using Test
using LinearAlgebra
using Random
using BSRandPivotQR

include(joinpath(@__DIR__, "..", "benchmark", "rand_test_matrix.jl"))
include("test_bsqr_rand.jl")
include("test_bsqr_rand_bounds.jl")
