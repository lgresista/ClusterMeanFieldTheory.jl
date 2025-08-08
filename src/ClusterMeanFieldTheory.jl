module ClusterMeanFieldTheory

using LatticeUtilities
using LinearAlgebra
using NLsolve
using SparseArrays
using KrylovKit

include("exact_diagonalization.jl")
include("mean_field_cluster.jl")
include("cluster_initialization.jl")

# exact diagonalization 
export HeisenbergInteraction, SpinCluster
export calculate_hamiltonianmatrix!, calculate_hamiltonianmatrix
export calculate_spinoperators
export expectation_value, crossproduct_expectation
export expectation_value_subsystem, get_subsystem_state_index, get_environment_groups
export eigenmin
export calculate_magnetizations!, calculate_magnetizations
export get_periodic_cluster_interactions, get_periodic_cluster

# CMFT
export MeanFieldCluster
export get_nsites, get_magnetic_fields
export set_magnetizations!, recalculate_magnetic_fields!
export calculate_groundstate_energy, calculate_energyshift, calculate_groundstate_and_energy
export fixedpoint_iteration!, anderson_acceleration!

# help cluster initialization
export UnitCell
export LabeledBond, get_bond, get_label
export get_meanfield_cluster_interactions, get_meanfield_cluster, calculate_nsites
export get_unitcell_3x, get_bonds_3x, loc_to_loc3x

end
