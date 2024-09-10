struct MeanFieldCluster
    spincluster :: SpinCluster
    meanfieldinteractions :: Vector{HeisenbergInteraction}
    magnetizations :: Vector{Vector{Float64}}    
end

# some getter functions for convenience
nsites(mfcluster :: MeanFieldCluster) = mfcluster.spincluster.nsites
magnetic_fields(mfcluster :: MeanFieldCluster) = mfcluster.spincluster.magnetic_fields
calculate_hamiltonianmatrix(mfcluster :: MeanFieldCluster) = calculate_hamiltonianmatrix(mfcluster.spincluster)
calculate_hamiltonianmatrix!(H, mfcluster :: MeanFieldCluster) = calculate_hamiltonianmatrix!(H, mfcluster.spincluster)

function set_magnetizations!(mfcluster :: MeanFieldCluster, new_magnetizations :: Vector{Vector{Float64}})
    for i in eachindex(mfcluster.magnetizations)
        mfcluster.magnetizations[i] .= new_magnetizations[i]
    end
    return nothing
end

# self-consistently converge mean-field cluster, iteratively updating the magnetic fields according to mean-field bonds
function converge!(mfcluster :: MeanFieldCluster; max_iterations = 100, abstol = 1e-8, reltol = 1e-4, verbose = true)
  
    if verbose println("Setting up self-consistent solution of meanfield cluster") end

    spinoperators = calculate_spinoperators(nsites(mfcluster))
    hamiltonianmatrix = calculate_hamiltonianmatrix(mfcluster.spincluster)
    groundstate_energy, groundstate = eigenmin(hamiltonianmatrix)
    new_magnetizations = calculate_magnetizations(spinoperators, groundstate)

    abserror = norm(new_magnetizations .- mfcluster.magnetizations)
    relerror = abserror/maximum(norm.(mfcluster.magnetizations))
    iteration = 0
    
    if verbose println("Starting iteration with initial absolute error = $(abserror) and relative error of $relerror") end

    while abserror > abstol && relerror > reltol && iteration < max_iterations

        # update magnetizations
        set_magnetizations!(mfcluster, new_magnetizations)

        # recalculate magnetic fields from new magnetizations
        recalculate_magnetic_fields!(mfcluster)

        # generate new Hamiltonian matrix with new fields
        calculate_hamiltonianmatrix!(hamiltonianmatrix, mfcluster.spincluster) #idea: update only parts coming from fields
        
        # diagonalize to obtain ground state
        groundstate .= eigenmin(hamiltonianmatrix)[2]

        calculate_magnetizations!(new_magnetizations, spinoperators, groundstate)
        
        abserror = norm(new_magnetizations .- mfcluster.magnetizations)
        relerror = abserror/maximum(norm.(mfcluster.magnetizations))
        iteration +=1
    end

    set_magnetizations!(mfcluster, new_magnetizations)
    recalculate_magnetic_fields!(mfcluster)

    is_converged = iteration < max_iterations
    if verbose println("Converged: $(is_converged). Iterations: $iteration/$(max_iterations). Absolute error: $abserror, relative error: $relerror") end

    return is_converged, iteration, abserror, relerror
end

# recalculate the effective magnetic fields from magnetizations and mean-field bonds
# NOTE: All real magnetic fields get set to zero. Not yet compatible!
function recalculate_magnetic_fields!(mfcluster :: MeanFieldCluster)

    # deref 
    fields = magnetic_fields(mfcluster)
    magnetizations = mfcluster.magnetizations
    
    # set to zero
    for field in fields
        field .= 0.0
    end
    
    #iterate over meanfield bonds to calculate new effective fields
    for interaction in mfcluster.meanfieldinteractions

        #unpack interaction
        i, j = interaction.sites
        J = interaction.J

        #add mean field field bond contributions J_ij * (<Si> S_j + S_i * <S_j> - <S_i><S_j>)
        fields[j] .+= J * magnetizations[i] # <S_i> S_j
        fields[i] .+= J * magnetizations[j] # S_i <S_j>
    end

    return nothing
end

# calculate constant contribution from mean-field bonds to the hamiltonian/energy
function calculate_energyshift(mfcluster :: MeanFieldCluster)
    # deref magnetizations
    magnetizations = mfcluster.magnetizations

    #initialize constant energy shift (due to mean-field bonds)
    energyshift = 0.0

    for interaction in mfcluster.meanfieldinteractions
        #unpack interaction
        i, j = interaction.sites
        J = interaction.J

        #add mean field field bond contributions -J_ij <S_i><S_j>  (important for energy comparisons)
        energyshift -= J * dot(magnetizations[i], magnetizations[j]) # <S_i><S_j>
    end

    return energyshift
end


# given a geometric unitcell, bonds, and the linear size of the spin cluster,
# calculate all inter and intracluster bonds. interclusterbonds will be treated
# in a mean-field fashion (i.e. added as magnetic fields), while intraclusterbonds
# will be treated exactly
function get_meanfield_cluster_interactions(uc, bonds, L)
    #initialize periodic lattice (containing all bonds)
    periodiclattice = Lattice(L = L, periodic = [true, true])

    #initialize open lattice (only containing intracluster bonds)
    openlattice = Lattice(L = L, periodic = [false, false])

    intercluster_interactions = HeisenbergInteraction[]
    intracluster_interactions = HeisenbergInteraction[]

    # iterate over all specified bonds
    for b in bonds
        
        bond = b[1]
        J = b[2]

        # get interactions out of neighbor tables
        
        # interactions insite the cluster
        intrabonds = filter(b -> !isempty(b), eachcol(build_neighbor_table(bond, uc, openlattice))) # filter in case table is empty
        
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
