module ClusterMeanFieldTheory

using LatticeUtilities
using LinearAlgebra
using NLsolve
using SparseArrays
using ArnoldiMethod

import LatticeUtilities: nsites

include("exact_diagonalization.jl")
include("mean_field_cluster.jl")

# exact diagonalization 
export HeisenbergInteraction, SpinCluster
export calculate_hamiltonianmatrix!, calculate_hamiltonianmatrix
export calculate_spinoperators
export expectation_value
export eigenmin
export calculate_magnetizations!,calculate_magnetizations

# CMFT
export MeanFieldCluster
export nsites, magnetic_fields
export set_magnetizations!, recalculate_magnetic_fields!
export calculate_groundstate_energy, calculate_energyshift
export get_meanfield_cluster_interactions
export fixedpoint_iteration!, anderson_acceleration!

end
