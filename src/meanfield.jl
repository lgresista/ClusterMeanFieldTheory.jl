function get_cluster_bonds(uc, bonds, L)
    L = [2, 2]

    periodiclattice = Lattice(L = L, periodic = [true, true])
    openlattice = Lattice(L = L, periodic = [false, false])

    locs = [site_to_loc(s, uc, periodiclattice) for s in 1:nsites(uc, periodiclattice)]
    pos = [loc_to_pos(l..., uc) for l in locs]

    interclusterbonds = []
    intraclusterbonds = []

    for b in bonds
        
        bond = b[1]
        J = b[2]

        intrabonds = eachcol(build_neighbor_table(bond, uc, openlattice))
        allbonds = eachcol(build_neighbor_table(bond, uc, periodiclattice))
        interbonds = setdiff(allbonds, intrabonds)
        
        for interbond in interbonds
            push!(interclusterbonds, (interbond, J))
        end
        
        for intrabond in intrabonds
            push!(intraclusterbonds, (intrabond, J))
        end
    end

    return interclusterbonds, intraclusterbonds
end

#function calculate_magnetic_fields(magnetizations, interclusterbonds)
#
#end