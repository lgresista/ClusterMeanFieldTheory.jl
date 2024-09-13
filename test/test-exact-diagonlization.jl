const σ = Matrix{ComplexF64}[[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
const id = [1. 0.; 0. 1.]

# construct spin operators as Kronecker product of Pauli matrices
function S(μ, i, N)
    ms = Matrix{ComplexF64}[id for _ in 1:N]
    ms[i] = σ[μ]/2
    return kron(ms...)
end

# get the magnetization/digits of the n-th basis state e_n in the Kronecker basis
function nkron_to_digits(n, N)
    digits = zeros(Int, N)
    for i in eachindex(digits)
        p = 2^(N-i+1)
        r = mod1(n, p)
        digits[i] = r <= p/2
    end 
    return digits
end

# do this for all basis states
function get_kron_digits(N)
    return [nkron_to_digits(n, N) for n in 1:2^N]
end

# map the Kronecker basis vectors to the ED basis vector (by a simple permutation)
toint(digits :: Vector{<:Int}) :: UInt = sum(digits[i] * 2^(i-1) for i in eachindex(digits))
get_kron_to_ed_permutation(N) = Int.(toint.(get_kron_digits(N))) .+ 1

## calculate Hamiltonian for 1D chain of spins with Heisenberg interactions and magnetic fields

# Kronecker basis version
function H_1D_chain_kron(Js, magnetic_fields)
    N = length(magnetic_fields)
    H = zeros(ComplexF64, 2^N, 2^N)
    
    spinoperators = [[S(μ, i, N) for μ in 1:3] for i in 1:N]

    for i in eachindex(Js)
        Si = spinoperators[i]
        Sj = spinoperators[i+1]
        H .+= sum(Js[i]*Si[μ]*Sj[μ] for μ in eachindex(Si))
    end
    
    for i in eachindex(magnetic_fields)
        Si = spinoperators[i]
        H .+= sum(magnetic_fields[i] .* Si)
    end
    return H
end

# ED Version
function H_1D_chain_ED(Js, magnetic_fields) # supply N-1 Heisenberg couplings and N magnetic fields
    positions = [[0.0, i] for i in eachindex(magnetic_fields)]
    interactions = [HeisenbergInteraction(Js[i], (i, i+1)) for i in eachindex(Js)]
    cluster = SpinCluster(positions, interactions, magnetic_fields)
    
    return calculate_hamiltonianmatrix(cluster)
end

@testset "Comparison to Kronecker product implementation" begin
    @testset "spin operators" begin
        for N in 2:9
            
            p = get_kron_to_ed_permutation(N)

            spinoperators_cmft = calculate_spinoperators(N)
            spinoperators_kron = [[Hermitian(S(μ, i, N)) for μ in 1:3] for i in 1:N]

            for Si in spinoperators_kron
                for μ in eachindex(Si)
                    Si[μ].data .= Si[μ].data[p, p]
                end
            end

            @test spinoperators_cmft ≈ spinoperators_kron    
        end
    end

    @testset "1D chain Hamiltonians" begin
        for N in 2:9
            for _ in 1:10
                N = 4
                Js = rand(N-1) .- 0.5
                magnetic_fields = [rand(3) .- 0.5 for _ in 1:N]
                p = get_kron_to_ed_permutation(N)

                @test H_1D_chain_ED(Js, magnetic_fields) ≈ H_1D_chain_kron(Js, magnetic_fields)[p, p]
            end
        end
    end
end
