using Printf, JuMP, DelimitedFiles, Mosek, MosekTools, LinearAlgebra

function read_data(INSTANCE_NAME, FORMULATION, DATA_PATH, REPO_PB)
  if FORMULATION == "NO_BLOCKS"
    BLOCKS=readdlm(joinpath(DATA_PATH,"decompositions","blocks_cholesky", INSTANCE_NAME*"_sdp_blocks.txt"))
    VAR = BLOCKS[:,2]
    B1 = ["B1" for i=1:size(BLOCKS,1)]
    BLOCKS = [B1 VAR]
    # CLIQUE_TREE = readdlm(joinpath(DATA_PATH, "cliquetree_MD",INSTANCE_NAME*"_sdp_cliquetree.txt"))
    CLIQUE_TREE = []
  else
    BLOCKS=readdlm(joinpath(DATA_PATH,"decompositions", "blocks_$FORMULATION",INSTANCE_NAME*"_sdp_blocks.txt"))
    CLIQUE_TREE = readdlm(joinpath(DATA_PATH,"decompositions", "cliquetree_$FORMULATION",INSTANCE_NAME*"_sdp_cliquetree.txt"))
  end
  CONSTRAINTS = readdlm(joinpath(DATA_PATH, REPO_PB, INSTANCE_NAME*"_sdp_ctr.txt"),skipstart=1)
  MATRICES = readdlm(joinpath(DATA_PATH, REPO_PB, INSTANCE_NAME*"_sdp_mat.txt"),skipstart=1)
  return BLOCKS, CLIQUE_TREE, CONSTRAINTS, MATRICES
end

#########################################################################################################################

function treat_constants(MATRICES_with_CONSTANTS, CONSTRAINTS)
  #find constant rows
  index_none_column2 = findall(x -> x == "NONE", MATRICES_with_CONSTANTS[:,2])
  index_none_column3 = findall(x -> x == "NONE", MATRICES_with_CONSTANTS[:,3])
  rows_to_remove = intersect(index_none_column2,index_none_column3)
  CONSTANTS = MATRICES_with_CONSTANTS[rows_to_remove,:]
  #remove constants from MATRICES
  MATRICES_without_CONSTANTS = MATRICES_with_CONSTANTS[setdiff(1:end, rows_to_remove), :]
  #define dict constant per objctr
  constants_per_objconstraint = Dict(CONSTRAINTS[i,1] => 0.0 for i=1:size(CONSTRAINTS, 1))
  constants_per_objconstraint["OBJ"] = 0.0
  for l=1:size(CONSTANTS,1)
    name = CONSTANTS[l,1]
    constants_per_objconstraint[name] = CONSTANTS[l,4]
  end
  return MATRICES_without_CONSTANTS, constants_per_objconstraint
end #function treat_constants

#########################################################################################################################
# liste des variables de chaque bloc
function create_list_of_variables_per_block(BLOCKS)
  DIFFERENT_BLOCKS = Set{String}(BLOCKS[:,1])
  block_var = Dict{String, Set{String}}()
  for block in DIFFERENT_BLOCKS
    block_var[block] = Set{String}()
  end
  for line=1:size(BLOCKS,1)
    block = BLOCKS[line,1]
    var = BLOCKS[line,2]
    block_var[block] = union(block_var[block],Set{String}([var]))
  end
  return block_var
end


#########################################################################################################################

function creating_and_solving_model(DIFFERENT_VAR, NB_BLOCKS, DIFFERENT_BLOCKS, block_var,CLIQUE_TREE, CONSTRAINTS, MATRICES, constants_per_constraint,my_solver, objscale)
  #initialize
  m = Model(with_optimizer(my_solver))
  NB_CONSTRAINTS = size(CONSTRAINTS,1)
  jumpX = Array{LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}}(undef,NB_BLOCKS)
  coeff_block = Dict{Tuple{String,String}, Set{String}}()
  mat_var = Dict{Tuple{String,String}, Dict{String,Any}}()
  println("Dictionary linking blocks with their set of variables created")
  for (block, var_list) in block_var
    for var1 in var_list
      for var2 in var_list
        mat_var[(var1,var2)] = Dict{String,Any}()
        coeff_block[(var1,var2)] = Set{String}()
      end
    end
  end
  println("Dictionary linking coeff with their couple of variables created")
  @printf("mat_var size is %8d\n", length(coeff_block))
  size_block = Dict{String, Int64}()
  i_block = 0
  for (block, var_list) in block_var
    i_block += 1
    i_var1 = 1
    size_block[block] = length(var_list)
    jumpX[i_block] =  @variable(m,[1:size_block[block],1:size_block[block]], base_name= "X_$block", PSD)
    for var1 in var_list
      i_var2 = 1
      for var2 in var_list
        mat_var[(var1,var2)][block] = jumpX[i_block][i_var1, i_var2]
        coeff_block[(var1,var2)] = union(coeff_block[(var1,var2)],Set{String}([block]))
        i_var2+=1
      end
      i_var1+=1
    end
  end
  println("Block variables defined in the model m : ", i_block, " - ", length(block_var))
  println("Dictionary linking variables, blocks and data coefficients done")
  xp = Dict{String,JuMP.GenericAffExpr}()
  xp["OBJ"] = 0*jumpX[1][1,1]
  for p = 1:NB_CONSTRAINTS
    xp[CONSTRAINTS[p,1]] = 0*jumpX[1][1,1]
  end
  NB_TERMS = size(MATRICES,1)
   #NO MEAN version
   for term = 1:NB_TERMS
     objctr = MATRICES[term,1]
     var1 = MATRICES[term,2]
     var2 = MATRICES[term,3]
     val = MATRICES[term,4]
     nb_block = length(coeff_block[(var1,var2)])
     min_size = Inf
     min_var = first(mat_var[(var1,var2)])[2]
     for (block, var) in mat_var[(var1,var2)]
       if length(block_var[block]) < min_size
         min_size = length(block_var[block])
         min_var = var
       end
     end
     add_to_expression!(xp[objctr], val*min_var)
   end
  println("Objective and constraints scalar products computed")
  my_timer = @elapsed @objective(m, Min , objscale*(xp["OBJ"]+constants_per_constraint["OBJ"]))
  @printf("%-35s%10.6f s\n", "Objective added to model m", my_timer)
  constraints_ref=Dict{String,JuMP.ConstraintRef}()
  for p = 1:NB_CONSTRAINTS
    name = CONSTRAINTS[p,1]
    if (CONSTRAINTS[p,2]=="EQ")
      constraints_ref[name] = @constraint(m,xp[name]+constants_per_constraint[name] == CONSTRAINTS[p,3])
    elseif (CONSTRAINTS[p,2]=="LEQ")
      constraints_ref[name] = @constraint(m,xp[name]+constants_per_constraint[name]<=CONSTRAINTS[p,4])
    elseif (CONSTRAINTS[p,2]=="GEQ")
      constraints_ref[name] = @constraint(m,xp[name]+constants_per_constraint[name]>=CONSTRAINTS[p,3])
    else
      constraints_ref[*(name,"_LEQ")] = @constraint(m,xp[name]+constants_per_constraint[name]<=CONSTRAINTS[p,4])
      constraints_ref[*(name,"_GEQ")] = @constraint(m,CONSTRAINTS[p,3]<=xp[name]+constants_per_constraint[name])
    end
  end
  println("initial constraints added to model m")
  #constraints linking common terms in blocks
  nb_coupling_constraints=0
  if NB_BLOCKS > 1
    for i in 1:size(CLIQUE_TREE,1)
      B1 = CLIQUE_TREE[i,1]
      B2 = CLIQUE_TREE[i,2]
      vars_B1 = block_var[B1]
      vars_B2 = block_var[B2]
      common_vars = [ var for var in intersect(vars_B1,vars_B2)]
      for i in 1:length(common_vars)
        var1 = common_vars[i]
        for j in i:length(common_vars)
          var2 = common_vars[j]
          JuMPvar1 = mat_var[(var1,var2)][B1]
          JuMPvar2 = mat_var[(var1,var2)][B2]
          constraints_ref[*("cc$(B1)_$(B2)_$var1","_$var2")]=@constraint(m,JuMPvar1-JuMPvar2==0)
          nb_coupling_constraints+=1
        end
      end
    end
  end
  println("nb_coupling_constraints primal : ", nb_coupling_constraints)
  println("nb variables primal : ", length(JuMP.all_variables(m)))
  for ctr_type in list_of_constraint_types(m)
    println("nb constraints type $ctr_type primal : ", num_constraints(m, ctr_type[1], ctr_type[2]))
  end
  # print(m)
  val, time, bytes, gctime, memallocs = @timed optimize!(m)
  # MOI.write_to_file(backend(m).optimizer.model, "dump.task.gz")
  mem = bytes*10^(-6)
  println("Resolution time : ", time)
  println("Memory consumption $mem Mo")
  println("\nObjective value: ", 1/objscale*JuMP.objective_value(m))
  Xoptimal=Array{Any}(nothing,NB_BLOCKS)
  for i=1:NB_BLOCKS
    Xoptimal[i]=JuMP.value.(jumpX[i])
  end
  return Xoptimal,JuMP.objective_value(m), JuMP.termination_status(m), time, mem
end

#########################################################################################################################

function sdp_model_solver_to_specify(INSTANCE_NAME, FORMULATION, DATA_PATH, REPO, my_solver, objscale)
  my_solver = Mosek.Optimizer
  #read data files
  my_timer=@elapsed (BLOCKS, CLIQUE_TREE, CONSTRAINTS, MATRICES_with_CONSTANTS)=read_data(INSTANCE_NAME, FORMULATION, DATA_PATH, REPO)
  @printf("%-35s%10.6f s\n", "read_data", my_timer)
  #treat constants
  my_timer = @elapsed (MATRICES, constants_per_objconstraint) = treat_constants(MATRICES_with_CONSTANTS, CONSTRAINTS)
  @printf("%-35s%10.6f s\n", "constant coeffcients treated", my_timer)

  #CONSTRAINTS
  DIFFERENT_CONSTRAINTS = Set{String}(CONSTRAINTS[:,1])
  NB_CONSTRAINTS = length(DIFFERENT_CONSTRAINTS)
  #variables
  DIFFERENT_VAR = Set{String}(BLOCKS[:,2])
  NB_VAR = length(DIFFERENT_VAR)
  #define block patterns
  DIFFERENT_BLOCKS = Set{String}(BLOCKS[:,1])
  NB_BLOCKS = length(DIFFERENT_BLOCKS)
  println("\n NB BLOCKS = $NB_BLOCKS \n")
  my_timer = @elapsed block_var = create_list_of_variables_per_block(BLOCKS)
  @printf("%-35s%10.6f s\n", "create_list_of_variables_per_block", my_timer)
  my_timer= @elapsed (Xoptimal,obj,status, time, mem) = creating_and_solving_model(DIFFERENT_VAR, NB_BLOCKS, DIFFERENT_BLOCKS, block_var,CLIQUE_TREE, CONSTRAINTS, MATRICES, constants_per_objconstraint,my_solver, objscale)
  #return S, Xoptimal, y
  @printf("%-35s%10.6f s\n", "creating_and_solving_model", my_timer)
  return Xoptimal, obj, status, time, mem
end

#########################################################################################################################

function solve_sdp(INSTANCE_NAME, FORMULATION, DATA_PATH)
    objscale = 4
    println("objscale : 10^(-$objscale) \n")
    Xoptimal, obj, status, time, mem = sdp_model_solver_to_specify(INSTANCE_NAME, FORMULATION, DATA_PATH,"OPF_instances", Mosek.Optimizer, 10.0^(-objscale))
    return
end
