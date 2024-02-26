using Combinatorics
using Gurobi
using JuMP
using CDDLib
using Polyhedra
using Random
using TickTock
using TimerOutputs




function main()

    # Create a TimerOutput, which we will use to keep track of 
    # the timing and memory usage of various lines of code
    Random.seed!(3)

    N_bar = 0
    E_bar = 0
    r = 0
    sales = 0
    
    # Create the output file
    f_out = open("../data/best_case_nested_R2.csv", "w")

    # Add the headers to the output file
    row_string = string("iter", ",",                      # Iteration that was run

                        "max_previous_assortment", ",",  # Expected revenue of best
														  # previously-offered assortment

                        "best_case_revenue", ",",         # Best-case expected revenue
                        "upper_bound",",",                # Upper bound
    					"n", ",",                         # Number of products
                        "time", ",",                      # Time
                        "rev_ordered"                     # true if past assortments are revenue ordered
						)              

    row_string = string(row_string, "\n");
    print(f_out, row_string);
    flush(f_out);

    for iter=1:100
       
        for setting ∈ ["rev_ordered","reverse_rev_ordered","other"]

            # Generate the past assortments
            if setting == "rev_ordered"
                # Revenue ordered assortments
                n = 10
                S = [[0,n]]
                for i=1:(n-1)
                    assortment = vcat([0],[k for k in (n-i):n])
                    push!(S,assortment)
                end

            elseif setting == "reverse_rev_ordered"
                # Reverse revenue ordered assortments
                n = 10
                S = [[0,n]]
                for i=1:(n-1)
                    assortment = vcat([k for k in 0:i],[n])
                    push!(S,assortment)
                end

            elseif setting == "other"
                # Complicated collection of past assortments
                S = [[0,3,8,13],[0,3,5,8,10,13,15],[0,1,3,5,6,8,10,11,13,15],[0,1,3,4,5,6,8,9,10,11,13,14,15],0:15]
            end    
    
            # Get the number of products and number of assortments
            M = length(S)
            n = maximum([maximum(S[m]) for m=1:M])
            K = 80
        
            # Get random revenues
            r_temp = 0
            while true
                r_temp = rand(1:10000,n)
                r_temp = sort(r_temp)
                if length(unique(r_temp)) == n
                    break
                end
            end
            r = Dict(i => r_temp[i] for i=1:n)
            r[0] = 0
    
            # Generate a random probability distribution
            Σ,λ_true = GetRandomDistribution(n,K)
   
            # Compute the historical sales data
            sales = GetFrequencies(λ_true,S,Σ)

            # Compute the best past assortment
            best_past_assortment = []
            best_past_assortment_obj_val = 0
            for m=1:M
                if sum(sales[m][i]*r[i] for i=S[m]) > best_past_assortment_obj_val
                    best_past_assortment = S[m]
                    best_past_assortment_obj_val = sum(sales[m][i]*r[i] for i=S[m])
                end
            end
 
            # Compute the optimal objective value for (RO)
            time_optimistic = @elapsed obj_val_best_case,optimal_assortment_best_case = SolveBestCase_Nested(S,sales,r)
               
            # Compute closed-form upper bound
            upper_bound = 0
            for k=1:K
    
                # Get the ranking
                σ = Σ[k]
    
                # Get the rank of the no-purchase option
                no_purchase_rank = σ[0]
    
                # Find the top ranked product
                for i=n:-1:1
                    if σ[i] < no_purchase_rank
                        upper_bound += r[i]*λ_true[k]
                        break
                    end
                end
            end
    
            @show obj_val_best_case, upper_bound
                    
            # Print to file
            row_string = string(iter, ",",
                                best_past_assortment_obj_val,",",
                                obj_val_best_case,",",
                                upper_bound, ",",
                                n, ",",
                                time_optimistic, ",",
                                setting)
    
            row_string = string(row_string, "\n");
            print(f_out, row_string);
            flush(f_out);
        end
    end
end      
main()
println("") 
