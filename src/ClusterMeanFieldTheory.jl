module ClusterMeanFieldTheory

using LatticeUtilities
using LinearAlgebra
using NLsolve
using SparseArrays
using KrylovKit
using Retry
using HDF5

import LatticeUtilities: nsites

include("exact_diagonalization.jl")
include("mean_field_cluster.jl")
include("cluster_initialization.jl")

# exact diagonalization 
export HeisenbergInteraction, SpinCluster
export calculate_hamiltonianmatrix!, calculate_hamiltonianmatrix
export calculate_spinoperators
export expectation_value
export eigenmin
export calculate_magnetizations!, calculate_magnetizations
export get_periodic_cluster_interactions, get_periodic_cluster

# CMFT
export MeanFieldCluster
export get_nsites, get_magnetic_fields
export set_magnetizations!, recalculate_magnetic_fields!
export calculate_groundstate_energy, calculate_energyshift
export fixedpoint_iteration!, anderson_acceleration!

# help cluster initialization
export LabeledBond, get_bond, get_label
export get_meanfield_cluster_interactions, get_meanfield_cluster, calculate_nsites
export get_unitcell_3x, get_bonds_3x, loc_to_loc3x
export run_cmft, run_batch_cmft

end
