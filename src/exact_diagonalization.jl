## Implement exact diagonalization for a spin model on a finite cluster with
## just Heisenberg interactions and unisotropic magnetic fields

# SU(2) symmetric Heisenberg Interaction
struct HeisenbergInteraction
    J::Float64
    sites::NTuple{2,Int}
end

# Struct capturing a spin model on a finite cluster 
struct SpinCluster
    positions::Vector{Vector{Float64}}
    interactions::Vector{HeisenbergInteraction}
    magnetic_fields::Vector{Vector{Float64}}
    nsites::UInt
end

# simplified constructor
function SpinCluster(
    positions::Vector{<:Vector{<:Real}},
    interactions::Vector{HeisenbergInteraction},
    magnetic_fields::Vector{<:Vector{<:Real}},
)
    return SpinCluster(positions, interactions, magnetic_fields, length(positions))
end

## functions to manipulate a state which is represented by an unsigned integer ##

# get bit at position i
function getbit(n::Unsigned, i::Integer)
    return (n & (one(n) << i)) != 0
end

# set bit at position i to 1
function setbit(n::Unsigned, i::Integer)
    return n | (one(n) << i)
end

# set bit at position i to 0
function unsetbit(n::Unsigned, i::Integer)
    return n & ~(one(n) << i)
end

# flip bit at position i
function flipbit(n::Unsigned, i::Integer)
    return n ⊻ (one(n) << i)
end

# flip bits as position i and j
function flipbits(n::Unsigned, i::Integer, j::Integer)
    n = flipbit(n, i)
    n = flipbit(n, j)
    return n
end

# convert bit (0/1) to spin (-0.5/0.5)
function getspin(bool)::Float64
    return -0.5 + bool
end

# get spin at position i
function getspin(n::Unsigned, i::Integer)::Float64
    return getspin(getbit(n, i))
end

## set up matrix representation of different operators ##

# calculate Hamiltonian matrix and store in H, given a spin model
function calculate_hamiltonianmatrix!(H::AbstractMatrix, spincluster::SpinCluster)
    N = spincluster.nsites

    H .= 0.0

    #iterate over all states
    for n in zero(UInt):(2^N - 1)
        # iterate over Heisenberg interactinos
        for interaction in spincluster.interactions
            # unpack interaction and shift to index begining at zero
            i = interaction.sites[1] - 1
            j = interaction.sites[2] - 1

            J = interaction.J

            if getbit(n, i) == getbit(n, j)
                #Sz contribution
                H[n + 1, n + 1] += J / 4
            else
                #Sz contribution
                H[n + 1, n + 1] -= J / 4

                #S± contribution (spin flip)
                m = flipbits(n, i, j)
                H[m + 1, n + 1] += J / 2
            end
        end

        #iterate over magnetic fields
        for i in 0:(N - 1)

            #get magnetic field vector at site i
            h = spincluster.magnetic_fields[i + 1]

            #get S^z eigenvalue (±1/2) of spin i
            si = getspin(n, i)

            #h^z S^z contribution
            H[n + 1, n + 1] += si * h[3]

            #h^x S^x + h^y S^y contribution
            m = flipbit(n, i)
            H[m + 1, n + 1] += h[1] / 2 + im * si * h[2]
        end
    end
end

# special version for Hermitian matrices
function calculate_hamiltonianmatrix!(H::Hermitian, spincluster::SpinCluster)
    return calculate_hamiltonianmatrix!(H.data, spincluster)
end

# out-of-place version (faster then initializing with spzeros for large matrices)
function calculate_hamiltonianmatrix(spincluster::SpinCluster)
    N = spincluster.nsites

    rows = UInt64[]
    cols = UInt64[]
    vals = ComplexF64[]

    #iterate over all states
    for n in zero(UInt):(2^N - 1)
        # iterate over Heisenberg interactinos
        for interaction in spincluster.interactions
            # unpack interaction and shift to index begining at zero
            i = interaction.sites[1] - 1
            j = interaction.sites[2] - 1

            J = interaction.J

            if getbit(n, i) == getbit(n, j)
                add_matrix_element!(rows, cols, vals, n + 1, n + 1, J / 4)
            else
                #Sz contribution
                add_matrix_element!(rows, cols, vals, n + 1, n + 1, -J / 4)

                #S± contribution (spin flip)
                m = flipbits(n, i, j)
                add_matrix_element!(rows, cols, vals, m + 1, n + 1, J / 2)
            end
        end

        #iterate over magnetic fields
        for i in 0:(N - 1)

            #get magnetic field vector at site i
            h = spincluster.magnetic_fields[i + 1]

            #get S^z eigenvalue (±1/2) of spin i
            si = getspin(n, i)

            #h^z S^z contribution
            add_matrix_element!(rows, cols, vals, n + 1, n + 1, si * h[3])

            #h^x S^x + h^y S^y contribution
            m = flipbit(n, i)
            add_matrix_element!(rows, cols, vals, m + 1, n + 1, h[1] / 2 + im * si * h[2])
        end
    end

    return Hermitian(sparse(rows, cols, vals, 2^N, 2^N))
end

# function to add entry to sparse matrix in coordinate format (COO)
function add_matrix_element!(rows, cols, vals, m, n, v)
    push!(rows, m)
    push!(cols, n)
    push!(vals, v)

    return nothing
end

function calculate_spinoperators(N)

    # spin operators as matrices for each site and each spin component S^x, S^y, S^z
    spinoperators = [
        Vector{Hermitian{ComplexF64,SparseMatrixCSC{ComplexF64,Int64}}}(undef, 3) for
        _ in 1:N
    ]

    # buffer for S^z eigenvalues
    szs = zeros(Float64, 2^N)

    # buffer state m from spinflip
    ms = zeros(UInt64, 2^N)

    #iterate over all sites
    for i in 0:(N - 1)

        #iterate over all states
        for n in zero(UInt):(2^N - 1)
            #get S^z eigenvalue (±1/2) of spin i
            szs[n + 1] = getspin(n, i)
            ms[n + 1] = flipbit(n, i) + 1
        end

        spinoperators[i + 1][1] = Hermitian(sparse(ms, 1:(2^N), fill(0.5, 2^N))) #x
        spinoperators[i + 1][2] = Hermitian(sparse(ms, 1:(2^N), im .* szs)) #y
        spinoperators[i + 1][3] = Hermitian(sparse(1:(2^N), 1:(2^N), szs)) #z
    end

    return spinoperators
end

## Diagonalization routines and observables ##

# calculate expectation value given a state and operator in matrix form
function expectation_value(operator::Hermitian, state)
    return real(dot(state, operator, state))
end

# calculate <state| S1 × S2 |state> for three-component operators S1 and S2
function crossproduct_expectation(S1, S2, state)
    v1 = S1 .* Ref(state)
    v2 = S2 .* Ref(state)

    return [
        dot(v1[2], v2[3]) - dot(v1[3], v2[2]),
        dot(v1[3], v2[1]) - dot(v1[1], v2[3]),
        dot(v1[1], v2[2]) - dot(v1[2], v2[1]),
    ]
end

# calculate the expectation value of an operator defined in a subsystem, i.e. of the form
# O_full = I_env ⊗ O_subsystem, by specifying O_subsystem and the sites of the subsystem.
function expectation_value_subsystem(
    operator::Hermitian,
    state,
    subsystem_sites,
    N;
    environment_groups=get_environment_groups(subsystem_sites, N),
)
    val = 0.0

    for n_env in eachindex(environment_groups)

        # basis states with fixed environment state (n_env - 1)
        state_indices = environment_groups[n_env]

        # project state into subsystem basis
        projected_state = state[state_indices]

        # calculate expectation value in subsystem
        val += expectation_value(operator, projected_state)
    end

    return val
end

# calculate projected state index into a subsystem with less sites
function get_subsystem_state_index(n::Integer, subsystem_sites::Vector{<:Integer})
    n_proj = zero(UInt)

    for (i, site) in enumerate(subsystem_sites)
        if getbit(n, site - 1) == true
            n_proj = setbit(n_proj, i - 1)
        end
    end
    return n_proj
end

# for a system with N sites and a subsystem with given sites, e.g. [2, 4, 5, 6]
# group all states by the state of their environment (complement of the subsystem)
# also saves the projected state index in the subsystem
function get_environment_groups(subsystem_sites, N)
    # get complement of subsystem
    environment_sites = [i for i in 1:N if i ∉ subsystem_sites]

    # number of sites in subsystem and environment
    N_s = length(subsystem_sites)
    N_e = length(environment_sites)

    # initialize projection table grouped by environment (one for each environment state)
    environment_groups = [zeros(Int, 2^N_s) for _ in 1:(2^N_e)]

    # loop over all global state indices
    for n in zero(UInt):(2^N - 1)
        # get state-index subsystem and environment
        n_proj = get_subsystem_state_index(n, subsystem_sites)
        n_env = get_subsystem_state_index(n, environment_sites)

        environment_groups[n_env + 1][n_proj + 1] = n + 1 # 
    end
    return environment_groups
end

# calculate lowest eigenvalue and vector using Lanczos method from KrylovKit.jl
function eigenmin(H::AbstractMatrix)
    vals, vecs, info = eigsolve(H, 1, :SR)
    return vals[1], vecs[1]
end

#= Old version not working on sparse matrices (faster only for small matrices)
function eigenmin(mat :: AbstractMatrix)
    vals, vecs = eigen(mat, 1:1)
    return vals[1], vecs[:, 1]
end
=#

#= old version using Arnoldi method from ArnoldiMethod.jl (compatible with sparse matrices but sometimes unstable)
function eigenmin(H :: Hermitian)
    # calculate partial schur decomposition
    decomp, history = partialschur(H, which = :SR, nev = min(size(H, 1), 9))

    # for Hermitian matrices the schur decomposition is an eigenvalue decomposition (otherwise add eigencomp(decomp))
    # select lowest eigenvalue (sometimes order not correct)
    idx = argmin(real.(decomp.eigenvalues))
    return real(decomp.eigenvalues[idx]), decomp.Q[:, idx]
end
=#

# calculate magnetization m_i^a = <state|S_i^a|state> for all spinoperators
# store into magnetizations
function calculate_magnetizations!(magnetizations, spinoperators, state)
    # iterate over sites
    for (m, S) in zip(magnetizations, spinoperators)
        # iterate over spin components (x, y, z)
        for i in eachindex(m)
            m[i] = expectation_value(S[i], state)
        end
    end
end

# out-of-place version
function calculate_magnetizations(spinoperators, state)
    magnetizations = [zeros(3) for _ in eachindex(spinoperators)]
    calculate_magnetizations!(magnetizations, spinoperators, state)
    return magnetizations
end