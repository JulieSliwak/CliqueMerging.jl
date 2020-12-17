include("..\\src\\solve_SDP.jl")

INSTANCE_NAME = "case1354pegase"
FORMULATION = "our_merging_heuristic"
DATA_PATH = "..\\data"

solve_sdp(INSTANCE_NAME, FORMULATION, DATA_PATH)
