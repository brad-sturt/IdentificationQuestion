using Combinatorics
using JuMP
using Gurobi
using Random

include("utils.jl")

function main(n=4)
  

	###################################################
 	# Initialize the output file, do other 
    # miscellaneous initializations
	###################################################
    
    # Set seed
    Random.seed!(1000)
    
    # Get the gurobi environment
    gurobi_env = Gurobi.Env()
    
    # Create the output file
    f_out = open("../data/appendix_A1.csv", "w")
    
    # Add the headers to the output file
    row_string = string("iter", ",",                      # Iteration that was run

                        "max_previous_assortments", ",",  # Expected revenue of best
														  # previously-offered assortment

    					"estimated_revenue", ",",         # Expected revenue under
														  # estimated ranking-based choice
														  # model of assortment found using
														  # estimate-then-optimize

    					"worst_case_revenue", ",",		  # Worst-case expected revenue

                        "best_case_revenue", ",",         # Best-case expected revenue

    					"n"                               # Number of products
						)              

    row_string = string(row_string, "\n");
    print(f_out, row_string);
    flush(f_out);


	###################################################
 	# Perform each iteration
	###################################################

	# Set the number of iterations 
    num_iterations = 1000

	# Iterate over the problem instances
    for iter=1:num_iterations
   

        ###################################################
     	# Initialize the problem instance
    	###################################################

        # Generate random revenues for the products
        revenues = rand(n)
        revenues = sort(revenues)
        r = Dict(i => revenues[i] for i=1:n)
        r[0] = 0
    
        # Construct the collection of revenue-ordered assortments
        rev_order = [v[1] for v in sort(collect(r), by = v -> v[2], rev=true)]   # rev_order[i] = product with i-th highest revenue, starting at 1
        revenue_ordered_assortments = [[0,rev_order[1]],]
        for i=2:n
            next_assortment = sort(vcat(revenue_ordered_assortments[i-1],[rev_order[i]]))
            push!(revenue_ordered_assortments,next_assortment)
        end
    
        # Generate a random `base' ranking-based choice model
        K = factorial(n+1)
        U = -1*log.(rand(factorial(n+1)))
        U_sum = sum(U)
        λ_base = U ./ U_sum

        # Construct the set of all distinct rankings Σ
        Σ_all = []
        for perm in permutations(0:n)
            push!(Σ_all, Dict(i-1 => perm[i] for i=1:(n+1)))
        end

        # Compute the historical sales data
        v = GetFrequencies(λ_base,revenue_ordered_assortments,Σ_all)


        ###################################################
     	# Evaluate the expected revenue under the best
        # previously-offered assortment
    	###################################################

        # Compute the expected revenue for each of the
        # past assortments
        expected_revenue_past_assortments = []
        for m=1:n
            rev = 0
            for i in revenue_ordered_assortments[m]
                rev += v[m][i]*r[i]
            end
            push!(expected_revenue_past_assortments, rev)
        end

        # Compute the expected revenue for the best past assortment
        max_expected_revenue_past_assortment,_ =  findmax(expected_revenue_past_assortments)
 

        ###################################################
     	# Find assortment using estimate-then-optimize
        ###################################################

        # Compute an estimate of the ranking-based choice model
        λ_est = EstimateDistribution(revenue_ordered_assortments,v,Σ_all,gurobi_env)

        # Compute the assortment which maximizes the expected revenue under the
        # estimated ranking-based choice model
        assortment_est, revenue_est = GetOptimalAssortment(λ_est,n,r,Σ_all,gurobi_env)
        assortment_est = sort(collect(assortment_est))


        ###################################################
     	# Evaluate the best-case and worst-case expected 
        # revenue of the assortment obtained using estimate-
        # then-optimize
        ###################################################
       
        revenue_wc,_ = EvaluateAssortment(assortment_est,r,revenue_ordered_assortments,v,n,Σ_all,gurobi_env,false)
        revenue_bc,_ = EvaluateAssortment(assortment_est,r,revenue_ordered_assortments,v,n,Σ_all,gurobi_env,true)
        println(revenue_est, ", ", assortment_est,", ", revenue_wc,", ", max_expected_revenue_past_assortment)


        row_string = string(iter, ",",
                            max_expected_revenue_past_assortment, ",",
        					revenue_est, ",",
        					revenue_wc, ",",
        					revenue_bc, ",",
        					n)
        row_string = string(row_string, "\n");
        print(f_out, row_string);
        flush(f_out);

    end
end
