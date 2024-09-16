# calculate mean-field cluster for the J1-J2 Heisenberg model on the square lattice.
# L specifies the size of the cluster (in the direction of each unit vector)
function get_mfcluster_square_j1j2_heisenberg(J, L)
    
    uc = UnitCell(lattice_vecs = [[1.,0.],[0.,1.]],
                         basis_vecs   = [[0.,0.]])

    bonds = [
        (Bond((1, 1), (1, 0)), J[1]),
        (Bond((1, 1), (0, 1)), J[1]),
        (Bond((1, 1), (1, 1)), J[2]),
        (Bond((1, 1), (1, -1)), J[2])
    ]
    
    
    intercluster_interactions, intracluster_interactions = get_meanfield_cluster_interactions(uc, bonds, L)
    
    lattice = Lattice(L = L, periodic = [true, true])
    locs = [site_to_loc(s, uc, lattice) for s in 1:nsites(uc, lattice)]
    pos = [loc_to_pos(l..., uc) for l in locs]
    spincluster = SpinCluster(collect.(pos), intracluster_interactions, [zeros(3) for _ in 1:length(pos)])
    mfcluster = MeanFieldCluster(spincluster, intercluster_interactions, [zeros(3) for _ in eachindex(pos)])    
    
    return mfcluster
end

# analytical expression for the ground-state energy of the dimer in the collinear phase
E0_dimer_collinear(J, s) = -J[1]/4 -J[1] * s^2 + 4 * J[2] * s^2 - sqrt(J[1]^2 + 4 * (J[1] - 4 * J[2])^2 * s^2)/2

@testset "Dimer cluster in collinear phase" begin
    for _ in 1:1000
        # random coupling
        J = (rand(2) .- 0.5) .* 3
        
        # initialize cluster
        dimer = get_mfcluster_square_j1j2_heisenberg(J, [2, 1])

        # random staggered magnetization
        s = rand() - 0.5
        
        # set exact magnetization and mean field magnetic fields
        m = normalize(rand(3)) .* s
        initial_magnetizations = [m, -m]
        set_magnetizations!(dimer, initial_magnetizations)
        recalculate_magnetic_fields!(dimer)
        
        # compare ground state energy with exact formula

        if !(calculate_groundstate_energy(dimer) ≈ E0_dimer_collinear(J, s)/2)
            @show J, s
        end

        @test calculate_groundstate_energy(dimer) ≈ E0_dimer_collinear(J, s)/2

    end
end