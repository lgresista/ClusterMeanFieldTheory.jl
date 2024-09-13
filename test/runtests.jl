using ClusterMeanFieldTheory
using ClusterMeanFieldTheory: setbit, unsetbit, getbit, flipbit, flipbits, getspin, HeisenbergInteraction, SpinCluster, calculate_hamiltonianmatrix, calculate_spinoperators
using LinearAlgebra
using Test

@testset "Bit-wise operations" begin
    include("test-bitwise-operations.jl")
end

@testset "Exact diagonalization" begin
    include("test-exact-diagonlization.jl")
end
