using ClusterMeanFieldTheory
using LinearAlgebra
using Test
using LatticeUtilities

import ClusterMeanFieldTheory as CMFT

@testset "Cluster Mean Field theory" begin
    @testset "Bit-wise operations" begin
        include("test-bitwise-operations.jl")
    end

    @testset "Exact diagonalization" begin
        include("test-exact-diagonlization.jl")
    end

    @testset "Mean field clusters of the J1-J2 Heisenberg model on the square lattice" begin
        include("test-square-lattice-clusters.jl")
    end

    @testset "Test 3x unit cell functions" begin
        include("test-uc-3x.jl")
    end
end