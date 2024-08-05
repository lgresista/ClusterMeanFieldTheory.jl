using ClusterMeanFieldTheory
using ClusterMeanFieldTheory: setbit, unsetbit, getbit, flipbit, flipbits, getspin
using Test

@testset "ClusterMeanFieldTheory.jl" begin
    @testset "Bit-wise-operations" begin
        include("test-bitwise-operations.jl")
    end
end
