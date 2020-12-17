include("..\\src\\generate_decompositions.jl")

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
