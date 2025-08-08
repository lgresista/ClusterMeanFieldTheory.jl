# ClusterMeanFieldTheory.jl

**ClusterMeanFieldTheory.jl** is a Julia package for solving quantum spin models using the cluster mean-field approximation (CMFT). In this approach, the system is divided into identical clusters: intra-cluster interactions are treated exactly via exact diagonalization, while inter-cluster couplings are approximated using a self-consistent mean-field scheme. For a detailed description of the method please refer to [Yong-Zhi Ren et al 2014 J. Phys.: Condens. Matter 26 115601](www.doi.org/10.1088/0953-8984/26/11/115601)

The package currently supports SU(2)-symmetric Heisenberg models on arbitrary spin clusters.

To install it simply write

```julia
]add https://github.com/lgresista/ClusterMeanFieldTheory.jl.git
```
In a julia terminal.

## Usage

In principle, one can define an arbitrary spin cluster by manually specifying positions of the cluster sites, interactions that are treated exactly and interaction that are to be mean-field approximated. For spin clusters with periodic boundary conditions that are multiples of the unit cell of conventional lattice spin model, however, we have automated this. We illustrate how to use the package for this usecase by showing the example of the $J_1-J_2$ Heisenberg model on the square lattice. 


First we have do define the cluster

```julia
using ClusterMeanFieldTheory
using LinearAlgebra

# The geometric unit cell is specified by defining a UnitCell (type reexported from LatticeUtilities.jl)
uc = UnitCell(lattice_vecs = [[1.,0.],[0.,1.]],
                     basis_vecs   = [[0.,0.]])

# implement nearest-neighbor interaction with label 1 
# and next-nearest neighbor interaction with label 2
bonds = [
    LabeledBond((1, 1), (1, 0), 1),
    LabeledBond((1, 1), (0, 1), 1),
    LabeledBond((1, 1), (1, 1), 2),
    LabeledBond((1, 1), (1, -1), 2)
]

# Define the cluster size
L = [2, 2]

# Define the coupling J = [J_1, J_2]
# first index is applied to bonds with label 1, second to bonds with label 2
J = [1.0, 0.2] 

# Obtain the mean-field cluster 
mfcluster = get_meanfield_cluster(uc, bonds, J, L);
```

The object ```mfcluster``` is of type ```MeanFieldCluster``` which has the following fields:
- ```mfcluster.magnetizations```: Current magnetizations $m_\alpha$
- ```mfcluster.meanfieldinteractions```: Interactions that are mean-field approximated
- ```mfcluster.spincluster```: All information for the single-cluster Hamiltonian

The object ```spincluster = mfcluster.spincluster``` is of type ```SpinCluster``` which saves the single-cluster Hamiltoniain $H_C$ that is used as input in the exact diagonalization. It has the follwing fields

- ```spincluster.positions```: Contains positions of the sites
- ```spincluster.interactions```: Contains the interactions of sites within the cluster that are treated exactly
- ```spincluster.magnetic_fields```: contains the effective fields $h_\alpha$
- ```spincluster.nsites```: Number of sites in the cluster

The function ```get_meanfield_cluster(uc, bonds, J, L)``` automatically determines all interactions within a cluster of size $L = L_1 \times L_2$ and saves them with the correct coupling $J_\mathrm{label}$ in the field ```mfcluster.spincluster.interactions```. It also determines the interactions that would occur between sites of different clusters (which will be mean-field approximated), and maps the sites back to sites within the cluster applying periodic boundary conditions. These are stored in ```mfcluster.meanfieldinteractions```. 

Next we choose initial magnetizations

```julia
# Choosing initial magnetizations (in this case randomly)
m0 = [rand(3).-0.5 for _ in 1:get_nsites(mfcluster)]

# Update these magnetizations in mfcluster
set_magnetizations!(mfcluster, m0)

# Recalculate the effective fields h_α (!!! not done automatically !!!)
recalculate_magnetic_fields!(mfcluster)
```
The function ```set_magnetizations!(mfcluster, m0)``` updates the initial magnetizations in ```mfcluster.magnetizations```, and ```recalculate_magnetic_fields!(mfcluster)``` recalculates the effective fields with the new magnetization and stores them in ```spincluster.magnetic_fields```.

Finally, we perform the fixed-point iteration
```julia
# Apply damped fixed points iterations
is_converged, iteration, abserror = fixedpoint_iteration!(
    mfcluster; # meanfield cluster
    max_iterations = 1000, # maximal number of iterations
    abstol = 1e-8, # tolerance after which iteration is stopped
    β = 0.5, # damping parameter
    verbose = false # toggle verbose output
)
```
The function ```fixedpoint_iteration!(mfcluster; max_iterations = 1000, abstol = 1e-8, β = 0.5, verbose = false)```  iteratively updates the magnetizations until they are converged. The converged magnetizations are stored in ```mfcluster.magnetizations```. 

Now, we can obtain the final ground state and its energy, and calculate other observables, as e.g. the spin-spin correlations:
```julia

# Final magnetization and ground state 
final_magnetizations = deepcopy(mfcluster.magnetizations)
e0, groundstate = calculate_groundstate_and_energy(mfcluster);

## Example for observables: Calculate correlations <ψ| S_i ⋅ S_j |ψ>

# Obtain matrix representation of all spin-operators in the format 
# [[S_1^x, S_1^y, S_1^z], [S_1^x, S_1^y, S_1^z]...]
spin_operators = calculate_spinoperators(mfcluster.spincluster.nsites)

# Select sites
i, j = 1, 2

Si = spin_operators[i]
Sj = spin_operators[j]

# Calculate <ψ| S_i ⋅ S_j |ψ> 
SiSj = sum(real(dot(groundstate, Si[μ], Sj[μ] * groundstate)) for μ in 1:3);
```

To obtain the ground-state energy and the ground state itself one can use ```calculate_groundstate_and_energy(mfcluster)```. The ground state is given in the local basis of the $S^z$ operator, where each site can have spin $s_i = (+1/2, -1/2) = (\uparrow, \downarrow) = (1, 0)$. The n-th basis state $|n\rangle = |s_1 s_2 s_3, \dots, s_N \rangle$ is then defined by the binary represenation $n = \sum_{i = 1}^N 2^{i-1} s_i$. 


For a general operator $\hat{O}$, expectation values can be calculated by ```dot(groundstate, O, groundstate)```, where $O$ is the matrix representation of the operator in the above basis. The matrix-representation of the spin-operators on each site and in each components for a cluster of site N can be obtained by  ```calculate_spinoperators(N)```. From these, the matrix representation of most operators can be calculated (although this may be costly for large clusters).
