# This script builds on the script "optimalPLD_exec_scripts/Toubia_Comparison.jl"
# from the Github repository associated with the following paper:
#
# D. Bertsimas and V. V. Mišić (2019). Exact first-choice product line 
# optimization. Operations Research, 67(3) : 651-670.

using CSV
using DataFrames
using Random
using Gurobi
using JuMP

include("general_utils.jl")
include("utils.jl")
include("construct_L.jl")
gurobi_env = Gurobi.Env()

function main()

    # Read in data
    revenues_mat = DataFrame(CSV.File("../conjoint_data/revenues_mat.csv",header=0))
    lambda_mat = DataFrame(CSV.File("../conjoint_data/lambda_mat.csv",header=0))
    orderings_mat = DataFrame(CSV.File("../conjoint_data/orderings_mat.csv",header=0))
    utilities_mat = DataFrame(CSV.File("../conjoint_data/utilities_mat.csv",header=0))

    # Set seed
    Random.seed!(3)

    # Create the output file
    f_out = open("../data/conjoint_R2.csv", "w")

    # Add the headers to the output file
    row_string = string("iter", ",",                      # Iteration
                        "max_previous_assortment", ",",   # Expected revenue of best past assortment
                        "best_case_revenue", ",",         # Optimal objective value of (OO)
                        "worst_case_revenue", ",",        # Optimal objective value of (RO)
                        "optimal_obj_val", ",",           # Optimal objective value of true problem
    					"n", ",",                         # Number of products
                        "M", ",",                         # Number of past assortments
                        "time_optimistic", ",",           # Time for solving (OO)
                        "time_robust", ",",               # Time for solving (RO)
                        "best_case_pareto", ",",          # Best-case expected revenue for current pareto assortment
                        "worst_case_pareto",",",          # Worst-case expected revenue for current pareto assortment
                        "expected_revenue_pareto"         # True expected revenue for current pareto assortment
						)              

    row_string = string(row_string, "\n");
    print(f_out, row_string);
    flush(f_out);

    for iter=1:10

        ################################################
        # Generate the subset of products to use
        # in the iteration
        ################################################

        # Choose a subset of the products
        num_products = 3584
        K = 330
        n = 15

        # Iterate over subsets of products until we find one
        # in which no two products have the same revenue
        while true
            prod_universe = randperm(num_products)[1:n]
            prod_universe = sort(prod_universe);
            push!(prod_universe, num_products+1);   
            r = collect(revenues_mat[1,:][prod_universe[1:n]])
            push!(r,0)
            if length(unique(r)) == n+1
                break
            end
        end        

        # Compute the correct order of the products
        prod_universe = prod_universe[sortperm(r)]
        r = sort(r)
        
        # Reindex the products
        prod_universe_dict = Dict()
        r_dict = Dict()
        for i=0:n
            prod_universe_dict[i] = prod_universe[i+1]
            r_dict[i] = r[i+1]
        end
        r = r_dict
        
        # Get a reverse mapping for products
        prod_universe_dict_inverse = Dict()
        for i=0:n
            prod_universe_dict_inverse[prod_universe_dict[i]] = i
        end
        
        
        # Give equal weight to each ranking
        λ = 1/K *ones(K)
        
        # Construct the set of rankings
        Σ = []
        for k=1:K
            # my understanding is that orderings_mat[k,l] is the l-th most preferred product
            # by customers of ranking k. Inverting it is necessary to have it match
            # the notation of our paper
            p2r = sortperm(collect(orderings_mat[k,:]))
            # Hence, p2r[i] is the rank of product i
        
            # Extract just the ranks of the relevant products
            p2r = p2r[prod_universe]
        
            # Get a mapping from product to relative rank
            p2rr = sortperm(sortperm(p2r))
            push!(Σ, Dict( i => p2rr[i+1]-1 for i=0:n))
        end


        ################################################
        # Compute the optimal assortment
        ################################################

        optimal_assortment, optimal_obj_val = GetOptimalAssortment(λ,n,r,Σ,gurobi_env)
 
        ################################################
        # Evaluate (RO) and (OO) for different choices
        # of numbers of past assortments
        ################################################
  
        for M ∈ 3:5


            ################################################
            # Generate the past assortments
            ################################################

            # Generate n-1 random integers between 1 and 2^(M+1) - 1  
            random_integers = rand(1:(2^(M+1)-1),n-1)
       
            # Generate M random past assortments
            S = [[0] for m=1:M]
            for i=1:(n-1)
                
                # If the m-th most significant bit of random_integer
                # is 1, then add product i to the m-th past assortment
                for m=1:M
                    if random_integers[i] % 2 == 1
                        push!(S[m],i)
                    end
                    random_integers[i] >>= 1
                end
            end
    
            # Add n-th product to all past assortments
            for m=1:M
                push!(S[m],n)
            end

    
            ################################################
            # Generate the data for past assortments
            ################################################

            # Compute the historical sales data
            sales = GetFrequencies(λ,S,Σ)

            # Compute the best past assortment
            best_past_assortment = []
            best_past_assortment_obj_val = 0
            for m=1:M
                if sum(sales[m][i]*r[i] for i=S[m]) > best_past_assortment_obj_val
                    best_past_assortment = S[m]
                    best_past_assortment_obj_val = sum(sales[m][i]*r[i] for i=S[m])
                end
            end
 
            ################################################
            # Solve (RO) and (OO)
            ################################################

            # Solve (OO)
            time_optimistic = @elapsed oo_obj_val, oo_optimal_assortment, oo_S_hat = SolveOptimisticOptimizationProblem(S, [i for i=0:n], sales, r, gurobi_env)

            # Solve (RO)
            time_robust =     @elapsed ro_obj_val, ro_optimal_assortment, ro_S_hat = SolveRobustOptimizationProblem(S, [i for i=0:n], sales, r, gurobi_env)


            ################################################
            # Evalaute the best-case and worst-case
            # expected revenue for each assortment in
            # the collection S_hat generated by (RO)
            ################################################

            # The following function outputs dictionaries that map
            # each assortment in ro_S_hat to the best-case 
            # and worst-case expected revenue for that assortment
            assortment_to_ro_obj_val, assortment_to_oo_obj_val = EvaluateWorstAndBestCases(ro_S_hat, S, [i for i=0:n], sales, r, gurobi_env)

            # Sort the assortments by their worst-case expected revenue,
            # and break ties by sorting by best-case expected revenue
            ro_S_hat = sort(collect(ro_S_hat), by=assortment -> (assortment_to_ro_obj_val[sort(collect(assortment))],assortment_to_oo_obj_val[sort(collect(assortment))]), rev=true)

            # Find the pareto assortments from the set ro_S_hat
            # and print to file
            prev_oo_obj_val = -Inf
            for assortment in ro_S_hat
                if prev_oo_obj_val + 0.01 < assortment_to_oo_obj_val[sort(collect(assortment))]

                    # `assortment` is a pareto assortment

                    # Update the bookkeeper for the current pareto assortment
                    prev_oo_obj_val = assortment_to_oo_obj_val[sort(collect(assortment))] 

                    # Compute the expected revenue for the current pareto assortment under
                    # the true choice model
                    sales_pareto = GetFrequencies(λ,[assortment],Σ)
                    expected_revenue_pareto = sum(sales_pareto[1][i]*r[i] for i=assortment)


                    # Evaluate the actual expected revenue of that assortment
                    row_string = string(iter, ",",
                                        best_past_assortment_obj_val,",",
                                        oo_obj_val,",",
                                        ro_obj_val,",",
                                        optimal_obj_val, ",",
                                        n, ",",
                                        M, ",",
                                        time_optimistic, ",",
                                        time_robust, ",",
                                        assortment_to_oo_obj_val[sort(collect(assortment))], ",", 
                                        assortment_to_ro_obj_val[sort(collect(assortment))], ",",
                                        expected_revenue_pareto
                                        )
           
                    row_string = string(row_string, "\n");
                    print(f_out, row_string);
                    flush(f_out);
                end
            end
        end
    end
end      
main()
println("") 


