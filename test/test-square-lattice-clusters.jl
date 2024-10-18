# analytical expression for the ground-state energy of the dimer in the collinear phase
E0_dimer_collinear(J, s) = -J[1]/4 -J[1] * s^2 + 4 * J[2] * s^2 - sqrt(J[1]^2 + 4 * (J[1] - 4 * J[2])^2 * s^2)/2

@testset "Dimer cluster in collinear phase" begin

    uc = UnitCell(lattice_vecs = [[1.,0.],[0.,1.]], basis_vecs= [[0.,0.]])

    bonds = [
        LabeledBond((1, 1), (1, 0), 1),
        LabeledBond((1, 1), (0, 1), 1),
        LabeledBond((1, 1), (1, 1), 2),
        LabeledBond((1, 1), (1, -1), 2)
    ]

    for _ in 1:1000
        # random coupling
        J = (rand(2) .- 0.5) .* 3
        
        # initialize cluster
        dimer = get_meanfield_cluster(uc, bonds, J, [2, 1])

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