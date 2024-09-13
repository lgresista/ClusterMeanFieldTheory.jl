struct MeanFieldCluster
    spincluster::SpinCluster
    meanfieldinteractions::Vector{HeisenbergInteraction}
    magnetizations::Vector{Vector{Float64}}
end

# some getter functions for convenience
nsites(mfcluster::MeanFieldCluster) = mfcluster.spincluster.nsites
magnetic_fields(mfcluster::MeanFieldCluster) = mfcluster.spincluster.magnetic_fields

# expand methods available for SpinCluster to MeanFieldCluster
function calculate_hamiltonianmatrix(mfcluster::MeanFieldCluster)
    return calculate_hamiltonianmatrix(mfcluster.spincluster)
end

function calculate_hamiltonianmatrix!(H, mfcluster::MeanFieldCluster)
    return calculate_hamiltonianmatrix!(H, mfcluster.spincluster)
end

# set magnetizations in MeanFieldCluster (DOES NOT UPDATE MAGNETIC FIELDS AUTOMATICALLY)
function set_magnetizations!(
    mfcluster::MeanFieldCluster, new_magnetizations::AbstractVector{<:AbstractVector}
)
    for i in eachindex(mfcluster.magnetizations)
        mfcluster.magnetizations[i] .= new_magnetizations[i]
    end
    return nothing
end

# set magnetizations with damping (DOES NOT UPDATE MAGNETIC FIELDS AUTOMATICALLY)
function set_magnetizations!(
    mfcluster::MeanFieldCluster, new_magnetizations::AbstractVector{<:AbstractVector}, β :: Float64
)
    for i in eachindex(mfcluster.magnetizations)
        m_old = mfcluster.magnetizations[i]
        m_new = new_magnetizations[i]
        @. mfcluster.magnetizations[i] = (1-β) * m_old + β * m_new
    end
    return nothing
end

# same functions for a flat magnetization vector instead of a vector of vectors (for NLsolve)
function set_magnetizations!(
    mfcluster::MeanFieldCluster, new_magnetizations::Vector{Float64}
)
    return set_magnetizations!(
        mfcluster, eachcol(reshape(new_magnetizations, (3, nsites(mfcluster))))
    )
end

function set_magnetizations!(
    mfcluster::MeanFieldCluster, new_magnetizations::Vector{Float64}, β
)
    return set_magnetizations!(
        mfcluster, eachcol(reshape(new_magnetizations, (3, nsites(mfcluster)))), β
    )
end


# recalculate the effective magnetic fields from magnetizations and mean-field bonds
# NOTE: All physical magnetic fields get set to zero. Not yet compatible with actual magnetic fields!
function recalculate_magnetic_fields!(mfcluster::MeanFieldCluster)

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
function calculate_energyshift(mfcluster::MeanFieldCluster)
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

function calculate_groundstate_energy(mfcluster :: MeanFieldCluster)    
    h = calculate_hamiltonianmatrix(mfcluster.spincluster)
    e0 = eigenmin(h)[1]
    eshift = calculate_energyshift(mfcluster)
    return (e0 + eshift)/nsites(mfcluster)
end


# given a geometric unitcell, bonds, and the linear size of the spin cluster,
# calculate all inter and intracluster bonds. interclusterbonds will be treated
# in a mean-field fashion (i.e. added as magnetic fields), while intraclusterbonds
# will be treated exactly
function get_meanfield_cluster_interactions(uc, bonds, L)
    #initialize periodic lattice (containing all bonds)
    periodiclattice = Lattice(; L=L, periodic=[true, true])

    #initialize open lattice (only containing intracluster bonds)
    openlattice = Lattice(; L=L, periodic=[false, false])

    intercluster_interactions = HeisenbergInteraction[]
    intracluster_interactions = HeisenbergInteraction[]

    # iterate over all specified bonds
    for b in bonds
        bond = b[1]
        J = b[2]

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

# self-consistently converge mean-field cluster, iteratively updating the magnetic fields according to mean-field bonds
function fixedpoint_iteration!(
    mfcluster::MeanFieldCluster; max_iterations=1000, abstol=1e-8, β = 1.0, verbose=true
)
    verbose && println("Setting up self-consistent solution of meanfield cluster")

    spinoperators = calculate_spinoperators(nsites(mfcluster))
    hamiltonianmatrix = calculate_hamiltonianmatrix(mfcluster.spincluster)
    groundstate_energy, groundstate = eigenmin(hamiltonianmatrix)
    new_magnetizations = calculate_magnetizations(spinoperators, groundstate)

    abserror = norm(new_magnetizations .- mfcluster.magnetizations)
    iteration = 0

    verbose && println("Starting iteration with initial absolute error = $(abserror)")

    while abserror > abstol && iteration < max_iterations

        # update magnetizations
        set_magnetizations!(mfcluster, new_magnetizations, β)

        # recalculate magnetic fields from new magnetizations
        recalculate_magnetic_fields!(mfcluster)

        # generate new Hamiltonian matrix with new fields
        calculate_hamiltonianmatrix!(hamiltonianmatrix, mfcluster) #idea: update only parts coming from fields

        # diagonalize to obtain ground state
        groundstate .= eigenmin(hamiltonianmatrix)[2]

        # calculate magnetization of ground-state
        calculate_magnetizations!(new_magnetizations, spinoperators, groundstate)

        abserror = norm(new_magnetizations .- mfcluster.magnetizations)
        iteration += 1
        #@show abserror
    end

    set_magnetizations!(mfcluster, new_magnetizations)
    recalculate_magnetic_fields!(mfcluster)

    is_converged = iteration < max_iterations
    verbose && println(
        "Converged: $(is_converged). Iterations: $iteration/$(max_iterations). Absolute error: $abserror",
    )

    return is_converged, iteration, abserror
end

# self-consistently converge mean-field cluster, iteratively updating the magnetic fields according to mean-field bonds
function anderson_acceleration!(
    mfcluster::MeanFieldCluster;
    max_iterations=1000,
    abstol=1e-8,
    m=5,
    beta=1,
    verbose=true,
    store_trace=false,
    show_trace=false,
    extended_trace=false
)
    verbose && println(
        "Setting up self-consistent solution of meanfield cluster via Anderson acceleration",
    )

    spinoperators = calculate_spinoperators(nsites(mfcluster))
    hamiltonianmatrix = calculate_hamiltonianmatrix(mfcluster.spincluster)
    groundstate_energy, groundstate = eigenmin(hamiltonianmatrix)
    new_magnetizations = calculate_magnetizations(spinoperators, groundstate)

    function f!(F, x)

        # update magnetizations
        set_magnetizations!(mfcluster, x)

        # recalculate magnetic fields from new magnetizations
        recalculate_magnetic_fields!(mfcluster)

        # generate new Hamiltonian matrix with new fields
        calculate_hamiltonianmatrix!(hamiltonianmatrix, mfcluster) #idea: update only parts coming from fields

        # diagonalize to obtain ground state
        groundstate .= eigenmin(hamiltonianmatrix)[2]

        # calculate magnetizations of ground state
        calculate_magnetizations!(new_magnetizations, spinoperators, groundstate)

        # parse result into output vector
        return F .= reduce(vcat, new_magnetizations)
    end

    # reshape magnetizations into flat vector
    x_initial = reduce(vcat, new_magnetizations)

    # find fixpoint via anderson acceleration
    sol = fixedpoint(
        f!,
        x_initial;
        iterations=max_iterations,
        ftol=abstol,
        m=m,
        beta = beta,
        store_trace=store_trace,
        show_trace=show_trace,
        extended_trace=extended_trace,
    )

    # update cluster with solution
    set_magnetizations!(mfcluster, sol.zero)
    recalculate_magnetic_fields!(mfcluster)

    return sol
end
