using Combinatorics
using Gurobi
using JuMP
using CDDLib
using Polyhedra
using Random
using TickTock
using TimerOutputs

include("utils.jl")
gurobi_env = Gurobi.Env()


###########################################################
# Construct the sets V and E
###########################################################

function ConstructVE(S,r)

    M = length(S)
    
    # Create the set B[1],…,B[M]
    B = Vector{Vector{Int64}}()
    push!(B,S[1])
    for m=1:(M-1)
        push!(B,[setdiff(Set(S[m+1]),Set(S[m]))...])
    end

    # Initialize the sets V_dict and E_dict
    V_dict = Dict{Int,Set{Tuple}}()
    E_dict = Dict{Int,Set{Tuple}}()
    
    # Populate V_dict[1]
    V_dict[1] = Set{Tuple}()
    for i in S[1]
        for j in S[1]
        #for κ in union(0,r[i])
            push!(V_dict[1], (i,r[j]))
            #push!(V_dict[1], (i,κ))
        end
    end

    # Populate remaining sets
    for m=1:(M-1)
        # Populate E_dict[m]
        E_dict[m] = Set{Tuple}()
        for (i,κ) in V_dict[m]
            for i_prime in union(i, B[m+1])
                #for κ_prime in union(κ, r[i_prime], [r[j] for j=B[m+1] if r[j] ≤ κ])
                for κ_prime in union(κ, [r[j] for j=B[m+1]])
                    push!(E_dict[m], ((i,κ),(i_prime,κ_prime)))
                end
            end
        end

        # Populate V_dict[m+1]
        V_dict[m+1] = Set{Tuple}()
        for ((i,κ),(i_prime,κ_prime)) in E_dict[m]
            push!(V_dict[m+1], (i_prime,κ_prime))
        end
    end

    # Combine the sets into N and E
    V = Set{Tuple}()
    E = Set{Tuple}()
    for m=1:M
        for (i,κ) in V_dict[m]
            push!(V, (m,i,κ))
        end
    end
    for m=1:(M-1)
        for ((i,κ),(i_prime,κ_prime)) in E_dict[m]
            push!(E, ((m,i,κ),(m+1,i_prime,κ_prime)))
        end
    end

    # Return the sets
    return V, E
end



###########################################################
# Solves the robust optimization problem when the past
# assortments are nested
###########################################################

function Nested(S,v,r,warm_start=[],fixed_assortment=[])

    # Get the number of past assortments
    M = length(S)

    # Get the number of products
    n = maximum([maximum(S[m]) for m=1:M])

    # Verify that the set of past assortments is nested
    for m=1:(M-1)
        if length(setdiff(Set(S[m]),Set(S[m+1]))) != 0
            error("S is not nested assortments")
        end
    end

    # Verify that the revenues are sorted and distinct
    if r[0] != 0
        error("r[0] is not equal to zero")
    end
    for i=0:(n-1)
        if r[i] ≥ r[i+1]
            error("The revenues for the products are not sorted and unique")
        end
    end

    # Create the set B[1],…,B[M]
    B = Vector{Vector{Int64}}()
    push!(B,S[1])
    for m=1:(M-1)
        push!(B,[setdiff(Set(S[m+1]),Set(S[m]))...])
    end

    # Create list of possible values of κ
    K_values = Set(r[i] for i=0:n)

    # Create the set of edges and nodes
    V, E = ConstructVE(S,r)
    counter = 0
    for ((m,i,κ),(m_prime,i_prime,κ_prime)) ∈ E
        if i == i_prime && κ == κ_prime
            counter += 1
        end
    end

    # Construct the optimization problem
    model = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_silent(model)

    # Create binary decision variables
    @variable(model, x[i=0:n], Bin)
    fix(x[0], 1)

    # If warm_start != [], then give warm start
    # values for the binary decision variables
    if length(warm_start) > 1
        for i=1:n
            if i in warm_start
                set_start_value(x[i],1)
            else
                set_start_value(x[i],0)
            end
        end
    end

    # If fixed_assortment != [], then give fixed
    # values for the binary decision variables
    if length(fixed_assortment) > 1
        for i=1:n
            if i in fixed_assortment
                fix(x[i],1)
            else
                fix(x[i],0)
            end
        end
    end

    # Create continuous decision variables
    @variable(model, α[m=1:M,i=S[m]])
    @variable(model, β[m=1:M,i=0:n,κ=K_values; (m,i,κ) ∈ V])
    @variable(model, γ[m=1:M,i=0:n,κ=K_values; (m,i,κ) ∈ V])

    # Create constraints
    expr1 = Dict{Any,AffExpr}()
    expr2 = Dict{Any,AffExpr}()
    expr3 = Dict{Any,AffExpr}()
    expr4 = Dict{Any,AffExpr}()
    expr5 = Dict{Any,AffExpr}()
    expr6 = Dict{Any,AffExpr}()
    for (m,i,κ) in V
        expr1[m,i,κ] = AffExpr(0.0)
        expr2[m,i,κ] = AffExpr(0.0)
        expr3[m,i,κ] = AffExpr(0.0)
        expr4[m,i,κ] = AffExpr(0.0)
        expr5[m,i,κ] = AffExpr(0.0)
        expr6[m,i,κ] = AffExpr(0.0)

        expr1[m,i,κ] += α[m,i]
        if m ∈ 1:(M-1)
            expr2[m,i,κ] -= β[m,i,κ]
        end
        if m ∈ 2:M
            expr3[m,i,κ] -= γ[m,i,κ]
        end
        if m == M
            expr4[m,i,κ] += κ
        end
        if i ∈ B[m] && κ != r[i]
            expr5[m,i,κ] += x[i]
        end
        if i ∉ B[m]
            for j=B[m]
                if κ == r[j]
                    expr5[m,i,κ] += x[i]
                end
            end
        end
        for j=B[m]
            if κ == r[j]
                expr6[m,i,κ] += 1 - x[j]
            end
        end
        @constraint(model, expr1[m,i,κ] + expr2[m,i,κ] + expr3[m,i,κ] ≤ expr4[m,i,κ] + r[n]*expr5[m,i,κ] + r[n]*expr6[m,i,κ])
    end

    @constraint(model, [((m,i,κ),(m_prime,i_prime,κ_prime)) in E], β[m,i,κ] + γ[m_prime,i_prime,κ_prime] ≤ 0)


    # Add objective function corresponding to method
    @objective(model, Max, sum(sum(v[m][i]*α[m,i] for i ∈ S[m]) for m=1:M))
    
    # Solve optimization problem
    optimize!(model)

    # Extract the optimal solution
    optimal_assortment = Vector{Int64}()
    for i=0:n
        if value(x[i]) > 0.5
            push!(optimal_assortment,i)
        end
    end
    return objective_value(model), optimal_assortment
end


function main()

    # Create a TimerOutput, which we will use to keep track of 
    # the timing and memory usage of various lines of code
    Random.seed!(3)

    V = 0
    E = 0
    r = 0
    sales = 0
    
    # Create the output file
    f_out = open("../data/appendix_A3.csv", "w")

    # Add the headers to the output file
    row_string = string("iter", ",",  # Iteration that was run

                        "n", ",",     # Number of products

                        "K", ",",     # Number of rankings with nonzero support
                       
                        "M", ",",     # Number of past assortments

                        "speed", ",", # Computation time in seconds

                        "max_previous_assortments", ",", # Expected revenue of best
														 # previously-offered assortment

    					"worst_case_revenue"             # Worst-case expected revenue of assortment
                                                         # found using algorithm from Section 4.2.1

                        )

    row_string = string(row_string, "\n");
    print(f_out, row_string);
    flush(f_out);

    ###################################################
 	# Perform each iteration for each value of n and M
	###################################################

    # Set the number of iterations, choice of K, and choices of
    # the number of products n
    num_iterations = 10
    Ks = [80,]
    ns = 20:20

    for n in ns, iter=1:num_iterations, K in Ks, M in 2:n

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
       
        # Decide how many elements to have in the first M-1 past assortments
        size_of_past_assortments = sort(randperm(n-1)[1:(M-1)])
        
        # Create a random permutation of the products
        rand_perm_of_products = randperm(n)

        # Create the assortments
        S = Vector{Vector{Int64}}()
        for m=1:(M-1)
            assortment = Vector{Int64}()
            push!(assortment, 0)
            for i in rand_perm_of_products[1:size_of_past_assortments[m]]
                push!(assortment, i)
            end
            sort!(assortment)
            push!(S, assortment)
        end
        push!(S, 0:n)

        # The following is the same as v, but in the correct format
        sales = GetFrequencies(λ_true,S,Σ)
        
        time_elapsed = @elapsed begin
 
            # For comparison, we find the best past assortment
            best_past_assortment = []
            best_past_assortment_obj_val = 0
            for m=1:M
                #println(S[m], ": ", sum(sales[m][i]*r[i] for i=S[m]))
                if sum(sales[m][i]*r[i] for i=S[m]) > best_past_assortment_obj_val
                    best_past_assortment = S[m]
                    best_past_assortment_obj_val = sum(sales[m][i]*r[i] for i=S[m])
                end
            end
     
            # Solve the robust optimization problem and extract the time 
            println("-------------------------------")
            computation_time = @elapsed obj_val_worst_case,optimal_assortment_worst_case = Nested(S,sales,r)
            @show computation_time, n,M, obj_val_worst_case,optimal_assortment_worst_case
        end

        # Print to file
        row_string = string(iter, ",",
                            n, ",",
                            K, ",",
                            M, ",",
                            time_elapsed, ",",
                            best_past_assortment_obj_val, ",",
                            obj_val_worst_case
                            )
 
        row_string = string(row_string, "\n");
        print(f_out, row_string);
        flush(f_out);
    end
end      
 
 main()
 println("") 
