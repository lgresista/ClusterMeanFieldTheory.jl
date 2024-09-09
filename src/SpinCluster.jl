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
function getspin(bool) :: Float64
    return -0.5 + bool
end

# get spin at position i
function getspin(n::Unsigned, i::Integer) :: Float64
    return getspin(getbit(n, i))
end

## set up matrix representation of different operators ##

# calculate Hamiltonian matrix and store in H, given a spin model
function calculate_hamiltonianmatrix!(H, spincluster::SpinCluster)
    N = spincluster.nsites
   
    #iterate over all states
    for n in zero(UInt):(N^2 - 1)
        # iterate over Heisenberg interactinos
        for interaction in spincluster.interactions

            #unpack interaction and shift to index begining at zero
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
        for i in 0:N-1

            #get magnetic field vector at site i
            h = spincluster.magnetic_fields[i+1]

            #get S^z eigenvalue (±1/2) of spin i
            si = getspin(n, i)

            #h^z S^z contribution
            H[n + 1, n + 1] += si * h[3]

            #h^x S^x + h^y S^y contribution
            m = flipbit(n, i)
            H[m + 1, n + 1] += h[1]/2 + im * si * h[2]
        end
    end
end

calculate_hamiltonianmatrix!(H :: Hermitian, spincluster :: SpinCluster) = calculate_hamiltonianmatrix!(H.data, spincluster)

# out-of-place version
function calculate_hamiltonianmatrix(spincluster::SpinCluster)
    N = spincluster.nsites
    H = Hermitian(zeros(ComplexF64, N^2, N^2))
    calculate_hamiltonianmatrix!(H, spincluster)
    return H
end


# needs testing!
# calcualte matrix representation of all spin operators S_i^a for N sites,
# where i = 1:N, a = x,y,z
function calculate_spinoperators(N)
    # spin operators as matrices for each site and each spin component S^x, S^y, S^z
    spinoperators = [[Hermitian(zeros(ComplexF64, N^2, N^2)) for _ in 1:3] for _ in 1:N]

    #iterate over all states
    for n in zero(UInt):(N^2 - 1)

        #iterate over all sites
        for i in 0:N-1

            #get S^z eigenvalue (±1/2) of spin i
            si = getspin(n, i)

            #get state with flipped bit/spin at site i
            m = flipbit(n, i)

            #add S^x, S^y and S^z contributions
            S = spinoperators[i+1]
            S[1].data[m + 1, n + 1] += 1/2
            S[2].data[m + 1, n + 1] += im * si
            S[3].data[n+1,n+1] += si
        end
    end

    return spinoperators
end    

## Diagonalization routines and observables ##

# calculate expectation value given a state and operator in matrix form
function expectation_value(operator :: Hermitian, state)
    return real(dot(state, operator, state))
end

# calculate the minmal eigenvalue and associated eigenvector of a matrix
function eigenmin(mat)
    vals, vecs = eigen(mat, 1:1)
    return vals[1], vecs[:, 1]
end

# calculate magnetization m_i^a = <state|S_i^a|state> for all spinoperators
# store into magnetizations
function calculate_magnetizations!(magnetizations, spinoperators, state)
    for (m, S) in zip(magnetizations, spinoperators)
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
