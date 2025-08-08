# ClusterMeanFieldTheory.jl

**ClusterMeanFieldTheory.jl** is a Julia package for solving quantum spin models using the cluster mean-field approximation (CMFT). In this approach, the system is divided into identical clusters: intra-cluster interactions are treated exactly via exact diagonalization, while inter-cluster couplings are approximated using a self-consistent mean-field scheme.

The package currently supports SU(2)-symmetric Heisenberg models on arbitrary spin clusters.

To install it simply write

```julia
]add https://github.com/lgresista/ClusterMeanFieldTheory.jl.git
```
In a julia terminal.

## The cluster-mean-field approximation

This package is designed to use the cluster-mean-field approxmiation to find the ground state of quantum spin models. Currently it is limited to Heisenberg Hamiltonians of the form

$$
    H = \sum_{ij} J_{ij} \mathbf{S}_i \mathbf{S}_j \, .
$$

with SU(2) spin-operators $\mathbf{S}_i$ with $S = 1/2$.

The main idea of CMFT is to devidie the full crystal lattice into small clusters, and then treat the interactions of spins within a small cluster $c$ exactly, but approximate interactions of spins in different clusters with a standard mean-field decoupling

$$
    \mathbf{S}_i  \mathbf{S}_j \approx \langle \mathbf{S}_i \rangle  \mathbf{S}_j + \mathbf{S}_i  \langle \mathbf{S}_j \rangle - \langle \mathbf{S}_i \rangle  \langle \mathbf{S}_j \rangle \, .
$$

This makes it possible to rewrite the Hamiltonian as a sum of single-cluster Hamiltonians

$$
H = \sum_{i} H_c
$$

of the form $

$$
H_c = \sum_{i} J_{ij} \mathbf{S}_i \mathbf{S}_j + h_i^C

$$

where the effective fields $\mathbf{h}_i^c$ and the energy shift $C^c$ originate from the inter-cluster mean-field interactions and a priori depend on the magnetizations $\langle \mathbf{S}_j \rangle$ of sites also in neighboring clusters. For the clusters to fully decouple, we use periodic boundary conditions by assuming that the magnetization patterns repeat identically across all clusters. To this end, we split the site index into $i = (c, \alpha)$, where $c$ denotes the cluster and $\alpha = 1, \dots, N_c$ the site within this cluster. Periodic boundary conditions then imply 
$$
\langle \mathbf{S}_{c\alpha} \rangle = \langle \mathbf{S}_{c'\alpha}\rangle \equiv \mathbf{m}_\alpha \, ,
$$

to hold for all clusters $c, c'$. This defines the $N_C$ cluster-independent magnetizations $\mathbf{m}_\alpha$, from which the now cluster-independent effective fields $\mathbf{h}_\alpha$ and the energy shift $C$ can be calculated using the formulas

$$
\mathbf{h_\alpha} = 
\sum_{c' \neq c} \sum_{\beta = 1}^{N_c}
\left(J_{c\alpha, c'\beta} + J_{c'\beta, c\alpha}\right)\mathbf{m}_\beta
$$

and the constant energy shift \( C \) is

$$
C = 
\sum_{c' \neq c} \sum_{\alpha,\beta = 1}^{N_c}
\left(J_{c\alpha, c'\beta} + J_{c'\beta, c\alpha}\right)
\mathbf{m}_\alpha \mathbf{m}_\beta \, 
$$

where $c$ is now an arbitrary reference cluster.

Since the single-cluster Hamiltonian $H_c$ depends on the magnetizations $\mathbf{m}_\alpha$ and these in turn depend on the ground state of $H_c$, the magnetizations must be determined self-consistently. This is what is implemented in this package.

## Fixed-point iteration for the ground state

For the self-consistent determination of the ground state, this package uses a damped fixed-point iteration, where at each iteration step $n$, an updated set of magnetizations $\mathbf{m}_\alpha^{n+1}$ is computed based on the values from the previous step $\mathbf{m}_\alpha^n$ as follows:  

1. Calculate the effective fields $\mathbf{h}_\alpha$ from $\mathbf{m}_\alpha^n$ using the expression for the effective fields given above.  
2. Determine the ground state of the resulting single-cluster Hamiltonian $H_c$ using the Lanczos algorithm.  
3. Calculate the magnetizations  
   $$
   \mathbf{m}_\alpha^{\mathrm{new}} = \langle \mathbf{S}_\alpha \rangle
   $$
   in this ground state.  
4. Update the magnetizations according to  
   $$
   \mathbf{m}_\alpha^{n+1} = (1 - \lambda) \mathbf{m}_\alpha^{\mathrm{new}} + \lambda \mathbf{m}_\alpha^n \, ,
   $$
   where $\lambda \in (0, 1]$ is a damping parameter used to improve convergence, typically set to $\lambda = 0.5$.  
5. Stop the iteration if  
   $$
   \sum_\alpha |\mathbf{m}_\alpha^{n+1} - \mathbf{m}_\alpha^n| < abstol \, ,
   $$
   otherwise continue with the next step.  

Once convergence is reached, one can calculate the final ground state using the self-consistent magnetizations. From this ground state, various observables can then be straightforwardly computed.

## Usage

In principle, one can define an arbitrary spin cluster by manually specifying positions of the cluster sites, interactions that are treated exactly and interaction that are to be mean-field approximated. For spin clusters with periodic boundary conditions that are multiples of the unit cell of conventional lattice spin model, however, we have automated this. We illustrate how to use the package for this usecase by showing the example of the $J_1-J_2$ Heisenberg model on the square lattice. We first provide the code and then explain what happens in detail below.

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

# Choosing initial magnetizations (in this case randomly)
m0 = [rand(3).-0.5 for _ in 1:get_nsites(mfcluster)]

# Update these magnetizations in mfcluster
set_magnetizations!(mfcluster, m0)

# Recalculate the effective fields h_α (!!! not done automatically !!!)
recalculate_magnetic_fields!(mfcluster)

# Apply damped fixed points iterations
is_converged, iteration, abserror = fixedpoint_iteration!(
    mfcluster; # meanfield cluster
    max_iterations = 1000, # maximal number of iterations
    abstol = 1e-8, # tolerance after which iteration is stopped
    β = 0.5, # damping parameter
    verbose = false # toggle verbose output
)

# Final results:
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

#### The mfcluster object

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

#### Setting initial magnetizations
The function ```set_magnetizations!(mfcluster, m0)``` updates the initial magnetizations in ```mfcluster.magnetizations```, and ```recalculate_magnetic_fields!(mfcluster)``` recalculates the effective fields with the new magnetization and stores them in ```spincluster.magnetic_fields```.

#### Performing the fixed point iteration
The function ```fixedpoint_iteration!(mfcluster; max_iterations = 1000, abstol = 1e-8, β = 0.5, verbose = false)``` which iteratively updates the magnetizations until they are converged. The converged magnetizations are stored in ```mfcluster.magnetizations```. 

#### Calculating ground-state observables
To obtain the ground-state energy and the ground state itself one can use ```calculate_groundstate_and_energy(mfcluster)```. The ground state is given in the local basis of the $S^z$ operator, where each site can have spin $s_i = (+1/2, -1/2) = (\uparrow, \downarrow) = (1, 0)$. The n-th basis state $|n\rangle = |s_1 s_2 s_3, \dots, s_N \rangle$ is then defined by the binary represenation $n = \sum_{i = 1}^N 2^{i-1} s_i$. 

For a general operator $\hat{O}$, expectation values can be calculated by ```dot(groundstate, O, groundstate)```, where $O$ is the matrix representation of the operator in the above basis. The matrix-representation of the spin-operators on each site and in each components for a cluster of site N can be obtained by  ```calculate_spinoperators(N)```. From these, the matrix representation of most operators can be calculated (although this may be costly for large clusters).