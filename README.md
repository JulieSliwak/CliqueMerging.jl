# CliqueMerging

This package contains a clique merging heuristic for OPF problems.

J. Sliwak, E. D. Andersen, M. F. Anjos, L. Létocart et E. Traversi, "A clique merging algorithm to solve semidefinite relaxations of optimal power flow problems," IEEE Transactions on Power Systems, vol. 36, no. 2, p. 1641-1644, 2021

## Required packages
* LinearAlgebra
* Printf
* DelimitedFiles
* JuMP
* Mosek
* MosekTools



## Data
All datasets are in the zip `data.zip`.

### MATPOWER data with more than 1000 buses
There are two files per instance:
* A ".mat" file with all quadratic terms (and possibly a constant for the objective function)
* A ".ctr" file that indicates the type of each constraint along with its lhs and rhs
NB : OPF problems with thermal limits in current, linear objective function and aggregation of generators associated to the same bus
* the files ".m" useful to reconstruct decompositions

### Clique decompositions
There are 15 clique decompositions in this repertory (10 decompositions used for linear regression and 5 obtained with our clique merging strategy). One clique decomposition is defined by two files:
* A block file
* A clique tree file

### Log files
There are 15 different log files for each instance:
* The ones used for the linear regression (10 decompositions)
* The ones obtained with our decomposition

### Linear regression parameters
There are the linear regression parameters used to generate our clique decomposition (1 file per instance).

## Functions

### To generate clique decompositions
```julia
write_file_decompositions(data_path, instances, algo, repo_name, FUSION, fusion_type = "", size_max_IPmerging=0, nb_times=0, kmax=0)
```
The function `write_file_decompositions` generates clique decompositions for instances in `instances` in repository `..\\data\\decompositions` following algorithm `algo` ("MD" (Minimum Degree), "cholesky", "cholesky_real") and clique merging algorithm `fusion_type` ("IP" (Integer Programming), "Molzahn", "our_merging_heuristic") if `FUSION="YES"`. The parameters `size_max_IPmerging` and `nb_times` are useful for the clique merging "IP": `size_max_IPmerging` is the clique size that we do not want to exceed and `nb_times` is the number of times that we want to apply the clique merging algorithm. The parameter `kmax` is useful for the clique merging "our_merging_heuristic" and indicates the number of iterations after which the criterion is changed. The argument `data_path` points the repertory which contains OPF instances. The argument `repo_name` is used to name the obtained decomposition e.g. "cholesky_cliquemergingIP".

Example:
```julia
data_path = "..\\data\\OPF_instances"
instances = ["case1354pegase.m"]
algo = "cholesky" #OR "MD" or "cholesky_real"
repo_name = "test_kmax2"
FUSION = "YES"
fusion_type = "our_merging_heuristic"
kmax=2
size_max_IPmerging = 0
nb_times = 0
write_file_decompositions(data_path, instances, algo, repo_name, FUSION, fusion_type, size_max_IPmerging, nb_times, kmax)
```

### To compute linear regression parameters
The file `linear_regression.jl` generates linear regression parameters from 10 clique decompositions.

### To solve an SDP problem
```julia
function solve_sdp(INSTANCE_NAME, FORMULATION, DATA_PATH)
```
The function `solve_sdp` solves the SDP problem for instance `INSTANCE_NAME` with clique decompositions specified in `FORMULATION`. The argument `DATA_PATH` points the repertory which contains OPF instances and clique decompositions.

Example:
```julia
INSTANCE_NAME = "case1354pegase"
FORMULATION = "our_merging_heuristic"
DATA_PATH = "..\\data"
solve_sdp(INSTANCE_NAME, FORMULATION, DATA_PATH)
```
