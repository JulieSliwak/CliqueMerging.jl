using LightGraphs, MetaGraphs, Xpress,  LinearAlgebra, SparseArrays, SuiteSparse

function load_matpower(filename)
  instance_name = split(filename, '.')[1]

  touch(instance_name*".temp")
  f = open(filename)
  out = open(instance_name*".temp", "w")

  # removing all ';' at end of lines
  while !eof(f)
    line = readline(f)
    if length(line) > 0 && line[1] != '%' && line[1] != 'f'
      s = split(line, ";")
      println(out, s[1])
    end
  end
  close(f)
  close(out)

  data = DelimitedFiles.readdlm(instance_name*".temp")
  rm(instance_name*".temp")
  return data
end

function find_numarray(i_start, data)
  i_debut = i_start
  while !isa(data[i_debut, 1], Int)
    i_debut+=1
  end
  i_fin=i_debut
  while !isa(data[i_fin,1], SubString)
    i_fin += 1
  end
  i_debut, i_fin-1
end

checkfor(data, line_ind, name) = (data[line_ind, 1] == name) || error("Expected ", name, " at line ", line_ind, ", got ", data[line_ind,1], " instead.")


function read_network_graph(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)
    if index_bus != Dict( i => i for i=1:nb_bus)
        println("!!! Specific numerotation of buses !!! \n")
    end
    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    G = MetaGraph(nb_bus)
    for i in 1:nb_bus
        set_props!(G, i, Dict(:name => "node$i"))
    end

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          LightGraphs.add_edge!(G, orig, dest)
      end
    end
    return G
end


function read_sparsity_pattern(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)

    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    sp = spzeros(nb_bus,nb_bus)

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          sp[orig,dest] = 1
      end
    end

    sp_sym = sp + sp'

    diag = zeros(nb_bus)

    for i=1:nb_bus
        sum_col = sum(sp_sym[i,:])
        diag[i] = Int(sum_col + 1)
    end
    return sp_sym + Diagonal(diag)
    # return sp+sp'+nb_bus*sparse(I, nb_bus, nb_bus)
end


function read_sparsity_pattern_real(instance_path::String)
    data = load_matpower(instance_path)
    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)

    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    sp = spzeros(2*nb_bus,2*nb_bus)

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          node_orig_Re = orig
          node_orig_Im = nb_bus + orig
          node_dest_Re = dest
          node_dest_Im = nb_bus + dest               #bus i => node i for iRe, node nb_bus+i for iIm
          sp[node_orig_Re, node_dest_Re] = 1
          sp[node_orig_Re, node_dest_Im] = 1
          sp[node_orig_Im, node_dest_Re] = 1
          sp[node_orig_Im, node_dest_Im] = 1
      end
    end
    sp_sym = sp + sp'
    diag = zeros(2*nb_bus)
    for i=1:2*nb_bus
        sum_col = sum(sp_sym[i,:])
        diag[i] = Int(sum_col + 1)
    end
    return sp_sym + Diagonal(diag)
end

function chordal_ext_cholesky(sparsity_pattern)
    A = sparsity_pattern
    nb_edges_A = (nnz(A) - size(A,1))/2
    #computing cholesky factorisation of A
    F = cholesky(A) #NOTE: order = F.p
    # println(F.p)
    # computing L + LT
    L = sparse(F.L)
    nb_edges_L = nnz(L) - size(A,1)
    nb_added_edges = nb_edges_L - nb_edges_A
    SP = L + L'
    #inverting permutation to get chordal extension of sparsity_pattern
    H = SP[invperm(F.p), invperm(F.p)]
    return H, F.p, nb_added_edges
end

function construct_graph_from_matrix(L)
    n = size(L,1)
    H = MetaGraph(n)
    for i in 1:n
        set_props!(H, i, Dict(:name => "node$i"))
    end
    for i in 1:n
        for j in 1:i
            if L[i,j] != 0
                MetaGraphs.add_edge!(H,i,j)
            end
        end
    end
    return H
end

##########################################################################################################

function MD_chordal_ext(G)
    n = nv(G)
    G_α = copy(G)
    H = copy(G)
    order = []
    nb_added_edges = 0
     for i=1:n
         degrees = [length(neighbors(G_α,s)) for s=1:(n-i+1)]
         perm = sortperm(degrees)
         v = perm[1]
         push!(order, Meta.parse(get_prop(G_α, v,:name)[5:end]))
         # println("v = ",v)
         # println(get_prop(G_α, v,:name))
         Nv = neighbors(G_α,v)
         for w in Nv
             for u in Nv
                 if w < u
                     if !has_edge(G_α, w, u)
                         add_edge!(G_α, w, u)
                         nodew = Meta.parse(get_prop(G_α, w,:name)[5:end])
                         nodeu = Meta.parse(get_prop(G_α, u,:name)[5:end])
                         add_edge!(H, nodew, nodeu)
                         nb_added_edges += 1
                    end
                     # println("edge = ",get_prop(G_α, w,:name),"--",get_prop(G_α, u,:name))
                end
            end
        end
        rem_vertex!(G_α,v)
    end
    return H, order, nb_added_edges
end

##########################################################################################################
function cliques_max(G, order)
    n = nv(G)
    C0 = []
    C = C0
    cliques_dict = Dict{Int64,Array{Any,1}}()
    k = 0
    for i=1:n
        v = order[i]
        Nv = neighbors(G,v)
        # println(v)
        C = [Nv[j] for j=1:length(Nv) if Nv[j] ∉ order[1:i]]
        C = [C ; v]
        # println(C)
        if C ⊈ C0
            k +=1
            cliques_dict[k] = C
            C0 = C
        end
    end
    return cliques_dict
end

##############################################################################################################

function weighted_graph(cliques_dict)
    nb_cl = length(cliques_dict)
    G = MetaGraph(nb_cl)
    for i in 1:nb_cl
        for j in 1:i
            inter = length(intersect(cliques_dict[i], cliques_dict[j]))
            if inter > 0
                LightGraphs.add_edge!(G,i,j)
                set_prop!(G, i, j, :weight, inter)
            end
        end
    end
    return G
end

##########################################################################################################

function Prim_algo(G)
    A = MetaGraph(nv(G))
    for i = 1:nv(G)
        set_prop!(A,i, :name, "C$i")
    end
    V = Set([ v for v in LightGraphs.vertices(G)])
    Vprim = Set([1])
    while Vprim != V
        max_weight = 0
        edge_to_add = LightGraphs.Edge(0,0)
        node_to_add = 0
        for edge in LightGraphs.edges(G)
            u = src(edge)
            v = dst(edge)
            if (u ∈ Vprim && v ∉ Vprim)
                weight_uv = get_prop(G,LightGraphs.Edge(u,v), :weight)
                if weight_uv > max_weight
                    max_weight = weight_uv
                    edge_to_add = LightGraphs.Edge(u,v)
                    node_to_add = v
                end
            elseif (u ∉ Vprim && v ∈ Vprim)
                weight_uv = get_prop(G,LightGraphs.Edge(u,v), :weight)
                if weight_uv > max_weight
                    max_weight = weight_uv
                    edge_to_add = LightGraphs.Edge(u,v)
                    node_to_add = u
                end
            end
        end
        LightGraphs.add_edge!(A, edge_to_add)
        set_prop!(A, edge_to_add, :weight, max_weight)
        push!(Vprim, node_to_add)
    end
    return A
end

##########################################################################################################

function mergingIP(cliques_dict, clique_tree,size_max)
    m = Model(with_optimizer(Xpress.Optimizer))
    x = Dict((Meta.parse(get_prop(clique_tree,src(e), :name)[2:end]),Meta.parse(get_prop(clique_tree,dst(e), :name)[2:end])) => @variable(m, base_name = "x_$(Meta.parse(get_prop(clique_tree,src(e), :name)[2:end]))_$(Meta.parse(get_prop(clique_tree,dst(e), :name)[2:end]))", binary=true) for e in edges(clique_tree))
    edge = collect(edges(clique_tree))[1]
    exp = Dict(i => @expression(m, 0*x[(Meta.parse(get_prop(clique_tree,src(edge), :name)[2:end]),Meta.parse(get_prop(clique_tree,dst(edge), :name)[2:end]))]) for i in keys(cliques_dict))
    for e in edges(clique_tree)
        C1 = Meta.parse(get_prop(clique_tree,src(e), :name)[2:end])
        C2 = Meta.parse(get_prop(clique_tree,dst(e), :name)[2:end])
        #x[(C1,C2)] = @variable(m, "x_$(C1)_$(C2)", Bin)
        size_merge_clique = length(union(cliques_dict[C1], cliques_dict[C2]))
        @constraint(m, size_merge_clique*x[(C1,C2)] <= size_max)
        add_to_expression!(exp[C1], 1, x[(C1,C2)])
        add_to_expression!(exp[C2], 1, x[(C1,C2)])
    end

    for (Ci, ex_ctr) in exp
        @constraint(m, exp[Ci]<=1)
    end
    @objective(m, Max, sum(length(intersect(cliques_dict[Meta.parse(get_prop(clique_tree,src(e), :name)[2:end])],cliques_dict[Meta.parse(get_prop(clique_tree,dst(e), :name)[2:end])]))*(length(intersect(cliques_dict[Meta.parse(get_prop(clique_tree,src(e), :name)[2:end])],cliques_dict[Meta.parse(get_prop(clique_tree,dst(e), :name)[2:end])]))+1)/2*x[(Meta.parse(get_prop(clique_tree,src(e), :name)[2:end]),Meta.parse(get_prop(clique_tree,dst(e), :name)[2:end]))] for e in edges(clique_tree)))
    optimize!(m)
    new_cliques_dict = copy(cliques_dict)
    new_clique_tree = copy(clique_tree)
    for (cliques,var) in x
        if value.(x[cliques]) == 1
            list_nodes1 = cliques_dict[cliques[1]]
            list_nodes2 = cliques_dict[cliques[2]]
            list_nodes_merge = union(list_nodes1, list_nodes2)
            new_cliques_dict[cliques[1]] = list_nodes_merge
            delete!(new_cliques_dict, cliques[2])
            v_iter = MetaGraphs.filter_vertices(new_clique_tree, :name, "C$(cliques[1])")
            v = 0
            for i in v_iter
                v = i
            end
            w_iter = MetaGraphs.filter_vertices(new_clique_tree, :name, "C$(cliques[2])")
            w = 0
            for j in w_iter
                w = j
            end
            N = neighbors(new_clique_tree, w)
            for n in N
                if n != v
                    LightGraphs.add_edge!(new_clique_tree, v, n)
                    weight = get_prop(new_clique_tree, LightGraphs.Edge(w,n), :weight)
                    set_prop!(new_clique_tree, v, n, :weight, weight)
                    # println("add_edge", get_prop(new_clique_tree,v, :name), get_prop(new_clique_tree,n, :name))
                end
            end
            rem_vertex!(new_clique_tree, w)
        end
    end
    new_list_maximal_cliques_clean = Dict{Int64, Array{Any,1}}()
    compteur = 1
    Dict_new_numbers = Dict{Int64, Int64}()
    for (clique, list_nodes) in new_cliques_dict
            new_list_maximal_cliques_clean[compteur] = list_nodes
            Dict_new_numbers[clique] = compteur
            compteur += 1
    end
     for i in 1:nv(new_clique_tree)
         old_name = get_prop(new_clique_tree, i, :name)
         new_number = Dict_new_numbers[parse(Int64,old_name[2:end])]
         set_prop!(new_clique_tree, i, :name, "C$new_number")
     end
     return new_list_maximal_cliques_clean, new_clique_tree
end

##########################################################################################################

function Molzahn_merging(cliques_dict, clique_tree, number_max_of_cliques)
    new_cliques_dict = copy(cliques_dict)
    new_clique_tree = copy(clique_tree)
    costs_fusion = Dict( e => cost_fusion(e, clique_tree, cliques_dict) for e in edges(clique_tree))
    while length(new_cliques_dict) > number_max_of_cliques
        (min_cost, edge_to_merge) = findmin(costs_fusion)
        clique1 = Meta.parse(get_prop(new_clique_tree, src(edge_to_merge), :name)[2:end])
        clique2 = Meta.parse(get_prop(new_clique_tree, dst(edge_to_merge), :name)[2:end])
        # println("fusion cliques $clique1 $clique2")
        list_nodes1 = new_cliques_dict[clique1]
        list_nodes2 = new_cliques_dict[clique2]
        list_nodes_fusion = union(list_nodes1,list_nodes2)
        # println("merged nodes : $list_nodes_fusion ")
        new_cliques_dict[clique1] = list_nodes_fusion
        delete!(new_cliques_dict, clique2)
        v_iter = MetaGraphs.filter_vertices(new_clique_tree, :name, "C$(clique1)")
        v = 0
        for i in v_iter
            v = i
        end
        w_iter = MetaGraphs.filter_vertices(new_clique_tree, :name, "C$(clique2)")
        w = 0
        for j in w_iter
            w = j
        end
        # println(w)
        N = neighbors(new_clique_tree, w)
        for n in N
            if n != v
                LightGraphs.add_edge!(new_clique_tree, v, n)
                weight = get_prop(new_clique_tree, LightGraphs.Edge(w,n), :weight)
                set_prop!(new_clique_tree, v, n, :weight, weight)
                # println("add_edge", get_prop(new_clique_tree,v, :name), get_prop(new_clique_tree,n, :name))
            end
        end
        rem_vertex!(new_clique_tree, w)
        costs_fusion = Dict(e => cost_fusion(e, new_clique_tree, new_cliques_dict) for e in edges(new_clique_tree))
    end
    new_list_maximal_cliques_clean = Dict{Int64, Array{Any,1}}()
    compteur = 1
    for (clique, list_nodes) in new_cliques_dict
            new_list_maximal_cliques_clean[compteur] = list_nodes
            compteur += 1
    end
    return new_list_maximal_cliques_clean
end


function cost_fusion(e, clique_tree, cliques_dict)
    clique1 = Meta.parse(get_prop(clique_tree, src(e), :name)[2:end])
    clique2 = Meta.parse(get_prop(clique_tree, dst(e), :name)[2:end])
    list_nodes1 = cliques_dict[clique1]
    list_nodes2 = cliques_dict[clique2]
    nintersect = length(intersect(list_nodes1,list_nodes2))
    elimmaxcliques1 = length(list_nodes1);
    elimmaxcliques2 = length(list_nodes2);
    lnewmaxcliques = elimmaxcliques1 + elimmaxcliques2 - nintersect
    nvarafter = (lnewmaxcliques)*(2*lnewmaxcliques+1) - (elimmaxcliques1*(2*elimmaxcliques1+1)+elimmaxcliques2*(2*elimmaxcliques2+1))
    ocostbefore = (nintersect)*(2*nintersect+1)
    cost = nvarafter - ocostbefore
    return cost
end

##########################################################################################################
#OUR DECOMPOSITION
##########################################################################################################
function cost_IPiteration(i, array_edges, cliques_dict, alpha, beta, m, nlc_before_merge)
    clique1 = array_edges[i,1]
    clique2 = array_edges[i,2]
    list_nodes1 = cliques_dict[clique1]
    # println(list_nodes1)
    list_nodes2 = cliques_dict[clique2]
    # println(list_nodes2)
    size_intersect = length(intersect(list_nodes1,list_nodes2))
    size_clique1 = length(list_nodes1)
    size_clique2 = length(list_nodes2)
    size_fusion_clique = size_clique1 + size_clique2 - size_intersect
    # println("Size clique 1 : $size_clique1, Size clique 2 : $size_clique2")
    #q1before_merge = sum(length(list_nodes)^3 for (clique, list_nodes) in cliques_dict) NOTE: inutile d'utiliser la somme totale car c'est la même pour tout le monde
    #q1 = q1before_merge - size_clique1^3 - size_clique2^3 + size_fusion_clique^3
    nlc = nlc_before_merge - size_intersect*(2*size_intersect+1)
    # println("nlc : $nlc")
    q2 = (m+nlc)^3
    cost = alpha*(-size_clique1^3-size_clique2^3+size_fusion_clique^3) + beta*q2
    return cost
end

function cost_Molzahn(i, array_edges, cliques_dict)
    clique1 = array_edges[i,1]
    clique2 = array_edges[i,2]
    list_nodes1 = cliques_dict[clique1]
    list_nodes2 = cliques_dict[clique2]
    nintersect = length(intersect(list_nodes1,list_nodes2))
    elimmaxcliques1 = length(list_nodes1);
    elimmaxcliques2 = length(list_nodes2);
    lnewmaxcliques = elimmaxcliques1 + elimmaxcliques2 - nintersect
    nvarafter = (lnewmaxcliques)*(2*lnewmaxcliques+1) - (elimmaxcliques1*(2*elimmaxcliques1+1)+elimmaxcliques2*(2*elimmaxcliques2+1))
    ocostbefore = (nintersect)*(2*nintersect+1)
    cost = nvarafter - ocostbefore
    return cost
end

function convert_clique_tree_to_array(A)
    array_edges = Array{Int64}(undef, ne(A),2)
    i = 0
    for e in edges(A)
        i += 1
        clique1 = Meta.parse(get_prop(A, src(e), :name)[2:end])
        clique2 = Meta.parse(get_prop(A, dst(e), :name)[2:end])
        array_edges[i,1] = clique1
        array_edges[i,2] = clique2
    end
    return array_edges
end

function fusion_Mosek_kmax_iterations(cliques_dict, clique_tree_array, number_max_of_cliques, alpha, beta, m, nlc, kmax)
    new_cliques_dict = copy(cliques_dict)
    array_edges = copy(clique_tree_array)
    nlc_before_merge = nlc
    costs_fusion = [cost_IPiteration(i, array_edges, new_cliques_dict, alpha, beta, m, nlc_before_merge) for i in 1:size(array_edges,1)]
    # println(costs_fusion)
    println("begin while...")
    it=0
    while length(new_cliques_dict) > number_max_of_cliques
        it+=1
        (min_cost, edge_to_merge) = findmin(costs_fusion)
        # println("Min cost : ", min_cost)
        # println("Edge to merge : ", edge_to_merge)
        clique1 = array_edges[edge_to_merge,1]
        clique2 = array_edges[edge_to_merge,2]
        # println("fusion cliques $clique1 $clique2")
        list_nodes1 = new_cliques_dict[clique1]
        list_nodes2 = new_cliques_dict[clique2]
        list_nodes_fusion = union(list_nodes1,list_nodes2)
        nintersect = length(intersect(list_nodes1,list_nodes2))
        # println("merged nodes : $list_nodes_fusion ")
        new_cliques_dict[clique1] = list_nodes_fusion
        delete!(new_cliques_dict, clique2)
        #delete edge in the clique tree array
        array_edges = array_edges[setdiff(1:end, edge_to_merge), :]
        #rename C2 into C1
        for i in 1:size(array_edges,1)
            if array_edges[i,1] == clique2
                array_edges[i,1] = clique1
            elseif array_edges[i,2] == clique2
                array_edges[i,2] = clique1
            end
        end
        # println(array_edges)
        nlc_before_merge -= (nintersect)*(2*nintersect+1)
        if it <= kmax
        # println("calculation new costs...")
            costs_fusion = [cost_IPiteration(i, array_edges, new_cliques_dict, alpha, beta, m, nlc_before_merge) for i in 1:size(array_edges,1)]
        # println(costs_fusion)
        else
            costs_fusion = [cost_Molzahn(i, array_edges, new_cliques_dict) for i in 1:size(array_edges,1)]
        end
    end
    new_list_maximal_cliques_clean = Dict{Int64, Array{Any,1}}()
    compteur = 1
    for (clique, list_nodes) in new_cliques_dict
            new_list_maximal_cliques_clean[compteur] = list_nodes
            compteur += 1
    end
    return new_list_maximal_cliques_clean
end

##########################################################################################################
#WRITE DECOMPOSITIONS
##########################################################################################################
function get_nb_linking_constraints(cliquetree)
    nb = 0
    for e in LightGraphs.edges(cliquetree)
        weight = get_prop(cliquetree, e, :weight)
        nb += weight*(2*weight+1)
    end
    return nb
end

function write_file_decompositions(data_path, instances, algo, repo_name, FUSION, fusion_type = "", size_max_IPmerging=0, nb_times=0, kmax=0)
    for instance in instances
        println(instance)
        instance_path = joinpath(data_path, instance)
        if algo == "MD_real"
            G = read_network_graph_real_numbers(instance_path)
        elseif algo == "cholesky" || algo == "Metis"
            sp = read_sparsity_pattern(instance_path)
        elseif algo == "cholesky_real"
            sp = read_sparsity_pattern_real(instance_path)
        else
            G = read_network_graph(instance_path)
        end
        if algo == "MD"
            H, order, nb_added_edges = MD_chordal_ext(G)
        elseif algo == "cholesky" || algo == "cholesky_real"
            L, order, nb_added_edges = chordal_ext_cholesky(sp)
            H = construct_graph_from_matrix(L)
        else
            println("error : algo $algo does not exist")
            exit()
        end
        nb_nodes = 0
        if algo == "cholesky" || algo == "cholesky_real"
            nb_nodes = size(sp,1)
        else
            nb_nodes = nv(G)
        end
        cliques_dict = cliques_max(H,order)
        GC = weighted_graph(cliques_dict)
        A = Prim_algo(GC)

        clique_tree = A
        if FUSION == "YES"
            if fusion_type == "IP"
                for i in 1:nb_times
                    new_cliques_dict, new_clique_tree = mergingIP(cliques_dict, clique_tree, size_max_IPmerging)
                    cliques_dict = new_cliques_dict
                    clique_tree = new_clique_tree
                end
            elseif fusion_type == "Molzahn"
                nb_max = 0.1*nb_nodes
                new_cliques_dict = Molzahn_merging(cliques_dict,A, nb_max)
                GC2 = weighted_graph(new_cliques_dict)
                clique_tree = Prim_algo(GC2)
                cliques_dict = new_cliques_dict
            elseif fusion_type == "our_merging_heuristic"
                nb_max = 0.1*nb_nodes
                nlc = get_nb_linking_constraints(A)
                infos_data_linear_regression = readdlm("..\\data\\parameters_linear_regression\\data_objscale10-4_$(instance[1:end-2]).csv", ';')
                m = infos_data_linear_regression[2,2]
                println(m)
                alpha = infos_data_linear_regression[13,1]
                println(alpha)
                beta = infos_data_linear_regression[13,2]
                println(beta)
                array_edges = convert_clique_tree_to_array(A)
                new_cliques_dict = fusion_Mosek_kmax_iterations(cliques_dict,array_edges, nb_max, alpha, beta, m, nlc, kmax)
                GC2 = weighted_graph(new_cliques_dict)
                clique_tree = Prim_algo(GC2)
                cliques_dict = new_cliques_dict
            end
        end

        isdir("..\\data\\decompositions\\blocks_"*repo_name) || mkpath("..\\data\\decompositions\\blocks_"*repo_name)
        f = open("..\\data\\decompositions\\blocks_"*repo_name*"\\$(instance[1:end-2])_sdp_blocks.txt", "w")
        if algo == "MD_real" || algo == "cholesky_real"
            nb_bus = Int(nb_nodes[instance]/2)
            for (clique, nodes_list) in cliques_dict
                for node in nodes_list
                    if node <= nb_bus
                        @printf(f, "%10s %20s\n", "B$clique", "VOLT_$(node)_Re")
                    else
                        bus = node - nb_bus
                        @printf(f, "%10s %20s\n", "B$clique", "VOLT_$(bus)_Im")
                    end
                end
            end
        else
            for (clique, nodes_list) in cliques_dict
                for node in nodes_list
                    @printf(f, "%10s %20s\n", "B$clique", "VOLT_$(node)_Re")
                    @printf(f, "%10s %20s\n", "B$clique", "VOLT_$(node)_Im")
                end
            end
        end
        close(f)

        isdir("..\\data\\decompositions\\cliquetree_"*repo_name) || mkpath("..\\data\\decompositions\\cliquetree_"*repo_name)
        f = open("..\\data\\decompositions\\cliquetree_"*repo_name*"\\$(instance[1:end-2])_sdp_cliquetree.txt", "w")

        if FUSION == "YES"
            for edge in LightGraphs.edges(clique_tree)
                clique1 = get_prop(clique_tree, src(edge), :name)[2:end]
                clique2 = get_prop(clique_tree, dst(edge), :name)[2:end]
                @printf(f, "%10s %20s\n", "B$clique1", "B$clique2")
            end
        else
            for edge in LightGraphs.edges(A)
                clique1 = src(edge)
                clique2 = dst(edge)
                @printf(f, "%10s %20s\n", "B$clique1", "B$clique2")
            end
        end
        close(f)
    end


    return
end
