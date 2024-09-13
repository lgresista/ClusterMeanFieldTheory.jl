module ClusterMeanFieldTheory
    
    using LatticeUtilities
    import LatticeUtilities: nsites
    using LinearAlgebra
    using NLsolve
    using SparseArrays

    include("exact_diagonalization.jl")
    include("mean_field_cluster.jl")

    export UnitCell, Bond, Lattice, nsites, site_to_loc, loc_to_pos
end
