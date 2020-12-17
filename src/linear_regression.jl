formulations = ["objscale10-4_BLOCKS_cholesky_ref0", "objscale10-4_BLOCKS_MD_ref0",
"objscale10-4_BLOCKS_cholesky_real_ref0", "objscale10-4_BLOCKS_cholesky_fusion_Dan_nb_max_10%_ref0", "objscale10-4_BLOCKS_cholesky_fusion50_x6_ref0",
"objscale10-4_BLOCKS_cholesky_fusion50_x5_ref0","objscale10-4_BLOCKS_cholesky_fusion50_x4_ref0", "objscale10-4_BLOCKS_cholesky_fusion50_x3_ref0",
"objscale10-4_BLOCKS_cholesky_fusion50_x2_ref0","objscale10-4_BLOCKS_cholesky_fusion50_x1_ref0"]


function read_mosek_log_more_info(log_file)
    lines = readlines(log_file)
    i = 1
    while split(lines[i], ':')[1]!= "nb_coupling_constraints primal " && split(lines[i], ':')[1]!="nb_coupling_constraints "
        i+=1
    end
    nlc = parse(Int, split(lines[i], ':')[2])
    while split(lines[i], ':')[1] != "  Constraints            "
        i+=1
    end
    n_total_constraints = parse(Int,split(lines[i], ':')[2])
    m = n_total_constraints - nlc
    while split(lines[i], ':')[1] != "Optimizer terminated. Time" && i < length(lines)
        i+=1
    end
    time = split(lines[i], ':')[2][2:end]
    nb_iter = split(lines[i-1], ' ')[1]
    objective = lines[i+2][18:end]
    return time, nb_iter, objective, m, nlc
end


nb_formulations = length(formulations)

DATA_PATH = joinpath("..","data","decompositions")

recap = open("data_linear_regression\\recap_linear_regression_objscale10-4.csv", "w")
recap_pareto = open("data_linear_regression\\recap_error_pareto_objscale10-4.csv", "w")
write(recap, "Instance ; alpha ; beta ; constant ; error ; with LOG10 ; Instance ; alpha ; beta ; constant ; error \n")
write(recap_pareto, "alpha * log10( sum |Ci|^3) + (1-alpha) * log10((m+nlc)^3) \n
Instance ; 0.1 ; 0.2 ; 0.3 ; 0.4 ; 0.5 ; 0.6 ; 0.7 ; 0.8 ; 0.9 \n")

for instance in instances
    println(instance)
    i = 0
    t = zeros(nb_formulations,1)
    q1 = zeros(nb_formulations,1)
    q2 = zeros(nb_formulations,1)
    e = ones(nb_formulations,1)
    errors_pareto = zeros(9,1)
    f = open("data\\parameters_linear_regression\\data_objscale10-4_$instance.csv", "w")
    write(f, "FORMULATION ; m ; nlc ; sum |Ci|^3 ; time per iteration  ;  ; log10( sum |Ci|^3) ;  log10((m+nlc)^3) ; ; alpha 0.1 ; alpha 0.2 ; alpha 0.3\n")
    for formulation in formulations
        # println(formulation)
        i += 1
        log_file = joinpath(pwd(), "..", "data", "logs", repo, "$(instance)_$(formulation).log")
        time = nb_iter = objective =  m =  nlc = "#"
        time_per_iteration = "#"
        remove_row = false
        try
            time, nb_iter, objective, m, nlc = read_mosek_log_more_info(log_file)
            time_per_iteration = parse(Float64,time)/parse(Int64,nb_iter)
            t[i,1] = time_per_iteration
            q2[i,1] = (m+nlc)^3
        catch
            println("WARNING : pb read log")
            remove_row = true
            time = nb_iter = objective =  m =  nlc = 0
            time_per_iteration = 0
            # t = t[1:end-1, 1]
            # q2 = q2[1:end-1, 1]
            # q1 = q1[1:end-1, 1]
            # e = e[1:end-1, 1]
        end
        #sum |C_i|
        #open blocks file
        if formulation[1:8] == "objscale"
            blocks_formulation = formulation[14:end]
        else
            blocks_formulation = formulation
        end
        if formulation[end-3:end-1] == "ref"
            blocks_formulation = blocks_formulation[1:end-5]
        end

        if blocks_formulation[1:end-3]=="BLOCKS_cholesky_fusion50"
            blocks_formulation = "blocks_cholesky_fusion_size_max_50_$(blocks_formulation[end-1:end])"
        end
        BLOCKS = readdlm(joinpath(DATA_PATH, blocks_formulation ,instance*"_sdp_blocks.txt"))
        block_var = create_list_of_variables_per_block(BLOCKS)
        sum_size_cliques_cube = 0
        for (block, list_var) in block_var
            size_block = length(list_var)
            sum_size_cliques_cube += size_block^3
        end
        if remove_row == false
            q1[i,1] = sum_size_cliques_cube
        else
            t=t[setdiff(1:end, i), :]
            q1=q1[setdiff(1:end, i), :]
            q2=q2[setdiff(1:end, i), :]
            e=e[setdiff(1:end, i), :]
            i-=1
        end

        write(f, "$formulation ; $m ; $nlc ; $sum_size_cliques_cube ; $time_per_iteration  ; ;  $(log10(sum_size_cliques_cube));  $(log10((m+nlc)^3)); ;")
         for i in 1:9
             alpha = 0.1*i
             criterion = alpha * log10(sum_size_cliques_cube) + (1 - alpha)*log10((m+nlc)^3)
             write(f, " $(criterion/time_per_iteration) ;")
             errors_pareto[i] += 0.5*(time_per_iteration-criterion)^2
        end
        write(f, "\n")
    end
    A = [q1'*q1 q1'*q2 q1'*e ; q2'*q1 q2'*q2 q2'*e ; e'*q1 e'*q2 e'*e ]
    b = [ q1'*t ; q2'*t ; e'*t]
    x = A\b
    error =  0.5*(t-(x[1]*q1+x[2]*q2+x[3]*e))'*(t-(x[1]*q1+x[2]*q2+x[3]*e))
    println("alpha : $(x[1]) ; beta : $(x[2]) ; constant : $(x[3]) \n")
    println("error 0,5 ||t-(alpha q1 + beta q2 + c)||^2 : ", error)
    write(f, "\n alpha ; beta ; constant \n $(x[1]) ; $(x[2]) ; $(x[3]) \n")
    write(f, "error 0,5 ||t-(alpha q1 + beta q2 + c)||^2 :  ; $(error[1]) \n")
    write(recap, "$instance ; $(x[1]) ; $(x[2]) ; $(x[3]) ; $(error[1]) ; ")
    write(f, "\n WITH LOG10 \n")
    q1 = log10.(q1)
    q2 = log10.(q2)
    A = [q1'*q1 q1'*q2 q1'*e ; q2'*q1 q2'*q2 q2'*e ; e'*q1 e'*q2 e'*e ]
    b = [ q1'*t ; q2'*t ; e'*t]
    x = A\b
    error =  0.5*(t-(x[1]*q1+x[2]*q2+x[3]*e))'*(t-(x[1]*q1+x[2]*q2+x[3]*e))
    println("alpha : $(x[1]) ; beta : $(x[2]) ; constant : $(x[3]) \n")
    println("error 0,5 ||t-(alpha log10(q1) + beta log10(q2) + c)||^2 : ", error)
    write(f, "\n alpha ; beta ; constant \n $(x[1]) ; $(x[2]) ; $(x[3]) \n")
    write(f, "error 0,5 ||t-(alpha log10(q1) + beta log10(q2) + c)||^2 :  ; $(error[1]) \n")
    close(f)
    write(recap, " WITHLOG10 ; $instance ; $(x[1]) ; $(x[2]) ; $(x[3]) ; $(error[1]) \n ")
    write(recap_pareto, "$instance ; ")
    for i in 1:9
        write(recap_pareto, "$(errors_pareto[i]) ; ")
    end
    write(recap_pareto, "\n")
end

close(recap)
close(recap_pareto)
