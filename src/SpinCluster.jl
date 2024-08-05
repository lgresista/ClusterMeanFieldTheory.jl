
struct HeisenbergInteraction
    J::Float64
    sites::NTuple{2,Int}
end

struct SpinCluster
    positions::Vector{Vector{Float64}}
    interactions::Vector{HeisenbergInteraction}
    magnetic_fields::Vector{Vector{Float64}}
    N::UInt
end

function SpinCluster(
    positions::Vector{<:Vector{<:Real}},
    interactions::Vector{HeisenbergInteraction},
    magnetic_fields::Vector{<:Vector{<:Real}},
)
    return SpinCluster(positions, interactions, magnetic_fields, length(positions))
end

function hamiltonianMatrix(spinCluster::SpinCluster)
    N = spinCluster.N
    H = zeros(ComplexF64, N^2, N^2)

    #iterate over all states
    for n in zero(UInt):(N^2 - 1)
        # iterate over Heisenberg interactinos
        for interaction in spinCluster.interactions

            #unpack interaction
            i, j = interaction.sites
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
            h = spinCluster.magnetic_fields[i+1]

            #get S^z eigenvalue (±1/2) of spin i
            si = getspin(n, i)

            #h^z S^z contribution
            H[n + 1, n + 1] += si * h[3]

            #h^x S^x + h^y S^y contribution
            m = flipbit(n, i)
            H[m + 1, n + 1] += h[1]/2 + im * si * h[2]
        end
    end

    return H
end