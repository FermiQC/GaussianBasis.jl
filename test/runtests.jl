using Base: _to_linear_index
using GaussianBasis
using Molecules
using HDF5
using Test

@testset verbose = true "GaussianBasis" begin
    include("test_integrals.jl")
    include("test_gradients.jl")
end