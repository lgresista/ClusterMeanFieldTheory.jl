## Mean field clusters from primitive unit cell and bonds
struct LabeledBond{D}
    bond::Bond{D}
    label::Int
end

function Base.show(io::IO, lbond::LabeledBond{D}) where {D}
    return print(
        io,
        "LabeledBond{$(D)}(displacement=$(lbond.bond.displacement), orbitals=$(lbond.bond.orbitals), label=$(lbond.label))",
    )
end

function LabeledBond(orbitals, displacement, label::Int)
    return LabeledBond(Bond(orbitals, displacement), label)
end

#Return a meanfield cluster from a primitive unitcell repeated L[i] times in each bravais dimension i.
function get_meanfield_cluster(
    uc::UnitCell, bonds::Vector{<:LabeledBond}, J::AbstractVector{<:Real}, L
)

    # generate inter- and intra-cluster interactions
    intercluster_interactions, intracluster_interactions = get_meanfield_cluster_interactions(
        uc, bonds, J, L
    )

    # calculate positions of lattice points
    lattice = Lattice(; L=L, periodic=[true, true])
    locs = [site_to_loc(s, uc, lattice) for s in 1:get_nsites(uc, L)]
    pos = [loc_to_pos(l..., uc) for l in locs]

    # initialize spin-cluster with intracluster interactions and zero magnetic field
    spincluster = SpinCluster(
        collect.(pos), intracluster_interactions, [zeros(3) for _ in 1:length(pos)]
    )

    # set up meanfield cluster with intercluster interactions (treated by mean-field) and zero magnetic field and magnetizations
    mfcluster = MeanFieldCluster(
        spincluster, intercluster_interactions, [zeros(3) for _ in eachindex(pos)]
    )

    return mfcluster
end

# given a geometric unitcell, bonds, couplings, and the linear size of the spin cluster,
# calculate all inter and intracluster bonds. interclusterbonds will be treated
# in a mean-field fashion (i.e. added as magnetic fields), while intraclusterbonds
# will be treated exactly
function get_meanfield_cluster_interactions(
    uc::UnitCell, bonds::Vector{<:LabeledBond}, Js, L
)
    #initialize periodic lattice (containing all bonds)
    periodiclattice = Lattice(; L=L, periodic=[true, true])

    #initialize open lattice (only containing intracluster bonds)
    openlattice = Lattice(; L=L, periodic=[false, false])

    intercluster_interactions = HeisenbergInteraction[]
    intracluster_interactions = HeisenbergInteraction[]

    # iterate over all specified bonds
    for b in bonds
        bond = b.bond
        J = Js[b.label]

        # get interactions out of neighbor tables

        # interactions insite the cluster
        intrabonds = filter(
            b -> !isempty(b), eachcol(build_neighbor_table(bond, uc, openlattice))
        ) # filter in case table is empty

        # all interactions
        allbonds = eachcol(build_neighbor_table(bond, uc, periodiclattice))

        # interactions between cluster is all bonds except the intraclusterbonds
        interbonds = setdiff(allbonds, intrabonds)

        for interbond in interbonds
            push!(intercluster_interactions, HeisenbergInteraction(J, Tuple(interbond)))
        end

        for intrabond in intrabonds
            push!(intracluster_interactions, HeisenbergInteraction(J, Tuple(intrabond)))
        end
    end

    return intercluster_interactions, intracluster_interactions
end

function get_nsites(uc, L)
    return uc.n * prod(L)
end

## 3x Unitcell: Get unitcell made up of 3 primitive unit cell on the locations (0, 0), (0, 1), (1, 0) (in bravais index notation)
# triple the basis vectors
function get_3x_basisvecs(uc)
    return collect.(
        [
            uc.basis_vecs
            uc.basis_vecs .+ Ref(uc.lattice_vecs[:, 1])
            uc.basis_vecs .+ Ref(uc.lattice_vecs[:, 2])
        ]
    )
end

# new bravais vectors
function get_latticevecs_3x(uc)
    lattice_vecs = collect.(eachcol(uc.lattice_vecs))
    return [lattice_vecs[1] .+ lattice_vecs[2], 2 .* lattice_vecs[1] .- lattice_vecs[2]]
end

# full new unitcell
function get_uc_3x(uc)
    return UnitCell(get_latticevecs_3x(uc), get_3x_basisvecs(uc))
end

# convert bravais indices
function ns_to_ns3x(n)
    return [n[1] / 3 + 2 / 3 * n[2], n[1] / 3 - n[2] / 3]
end

# and back
function ns3x_to_ns(nx)
    return [nx[1] + 2 * nx[2], nx[1] - nx[2]]
end

# convert location in terms of uc, to location in terms of uc3x
function loc_to_loc3x(loc, uc)
    n, o = loc

    ## find out in which of the three copied unitcells the loc lies

    # define copied unitcells
    uc_ns = ([0.0, 0.0], [1.0, 0.0], [0.0, 1.0])

    # calculate bravais indices on uc3 lattice assuming it lies in either of the three copied unit cells
    nxs = [round.(ns_to_ns3x(n .- uc_n), digits=10) for uc_n in uc_ns]

    # check for which of the copied unit cells the bravais indices are actual integers
    # --> it is an actual bravais vector
    uc3_idx = findfirst([all(isinteger.(nx)) for nx in nxs])

    # select the corresponding bravais indices 
    nx = nxs[uc3_idx]

    # find the basis index depending on in which copy of the unit cell the loc lives
    ox = o + (uc3_idx - 1) * length(uc.basis_vecs)

    return (nx, ox)
end

# convert bond of primitive unit cell to 3x unit cell
function get_bonds_3x(bonds, uc)
    uc_ns = ([0.0, 0.0], [1.0, 0.0], [0.0, 1.0])

    bonds3x = LabeledBond[]

    for uc_n in uc_ns, b in bonds
        bond = b.bond
        label = b.label

        # get locations for original unit cell
        loc_i = (uc_n, bond.orbitals[1])
        loc_j = (loc_i[1] .+ bond.displacement, bond.orbitals[2])

        # get locations for 3x extended unit cell
        loc3x_i = loc_to_loc3x(loc_i, uc)
        loc3x_j = loc_to_loc3x(loc_j, uc)

        # calculate displacement and orbitals for extended unit cell
        l = loc3x_j[1] .- loc3x_i[1]
        o = (loc3x_i[2], loc3x_j[2])

        # initialize new bond and push to bond3x vector
        bond3x = Bond(; displacement=l, orbitals=o)
        push!(bonds3x, LabeledBond(bond3x, label))
    end

    return bonds3x
end