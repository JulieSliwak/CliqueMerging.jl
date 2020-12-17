# CliqueMerging

This package contains a clique merging heuristic for OPF problems.

## Required packages
* LinearAlgebra
* Printf
* DelimitedFiles
* JuMP
* Mosek
* MosekTools



## Data
### MATPOWER data with more than 1000 buses
Two files per instance:
* one file "mat" with all quadratic terms (and possibly a constant for the objective function)
* one file "ctr" that indicates the type of each constraint along with its lhs and rhs
NB : OPF problems with thermal limits in current, linear objective function and aggregation of generators associated to the same bus
+ the files ".m" useful to reconstruct decompositions


### Clique decompositions
* Blocks files
* Clique tree files
15 clique decompositions in this repertory


### Log files
* The ones used for the linear regression (10 decompositions)
* The ones obtained with our decomposition


### Linear regression parameters
* Final linear regression parameters used to generate our clique decomposition (1 file per instance)

## Functions
### To solve an SDP problem
```julia
function solve_sdp(INSTANCE_NAME, FORMULATION, DATA_PATH)
```
The function `solve_sdp` solves the SDP problem for instance `INSTANCE_NAME` with clique decompositions specified in `FORMULATION`. The argument `DATA_PATH` point the repertory which contains OPF instances and clique decompositions.

### To generate clique decompositions

### To find linear regression parameters

## Example
