module ClusterMeanFieldTheory
    
    using LatticeUtilities
    import LatticeUtilities: nsites
    using LinearAlgebra
    

    include("SpinCluster.jl")
    include("MeanfieldCluster.jl")

    export UnitCell, Bond, Lattice, nsites, site_to_loc, loc_to_pos
end
