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
# The following function can run in several different ways,
# based on the choice of the `method` paramater.
#  - method = 0 will find an assortment that maximizes the 
#               best-case expected revenue 
#  - method = 1 will find an assortment that maximizes the
#               worst-case expected revenue
#  - method = 2 will find an assortment that maximizes the
#               best-case expected revenue subject to
#               the constraint that the worst-case expected
#               revenue is greater than or equal to the
#               argument `lower_bound`
# Moreover, if the `fixed_assortment` argument is populated, then the assortment
# will be fixed (and thus this function will serve as an evaluation
# of either the best-case or worst-case expected revenue, depending on
# `method`)
###########################################################

function Pareto_Nested(S,v,r,method,lower_bound=0,warm_start=[],fixed_assortment=[])

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

    # Create continuous decision variables for the worst-case reformulation
    @variable(model, α[m=1:M,i=S[m]])
    @variable(model, β[m=1:M,i=0:n,κ=K_values; (m,i,κ) ∈ V])
    @variable(model, γ[m=1:M,i=0:n,κ=K_values; (m,i,κ) ∈ V])

    # Create constraints for the worst-case reformulation
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

    # Create the decision variables for best-case
    @variable(model, f[m=1:(M-1),i=0:n,κ=K_values,i_prime=0:n,κ_prime=K_values; ((m,i,κ),(m+1,i_prime,κ_prime)) in E] ≥ 0)
    @variable(model, g[m=1:M,i=0:n,κ=K_values; (m,i,κ) ∈ V] ≥ 0)

    # Create primal constraints
    @constraint(model, [m=1:M,i=S[m]], sum(g[m_temp,i_temp,κ] for (m_temp,i_temp,κ) in V if m_temp == m && i_temp == i) == v[m][i])
    @constraint(model, [(m,i,κ) ∈ V; m ∈ 1:(M-1)], sum(f[m,i,κ,i_prime,κ_prime] for ((m_prime_prime,i_prime_prime,κ_prime_prime),(m_prime,i_prime,κ_prime)) in E if m_prime_prime == m && i_prime_prime == i && κ_prime_prime == κ) == g[m,i,κ])
    @constraint(model, [(m,i,κ) ∈ V; m ∈ 2:M], sum(f[m-1,i_prime,κ_prime,i,κ] for ((m_prime,i_prime,κ_prime),(m_prime_prime,i_prime_prime,κ_prime_prime)) in E if m_prime_prime == m && i_prime_prime == i && κ_prime_prime == κ) == g[m,i,κ])

    for (m,i,κ) ∈ V
        if i ∈ B[m] && κ != r[i]
            @constraint(model, g[m,i,κ] ≤ 1 - x[i])
        end
        if i ∉ B[m]
            for j=B[m]
                if r[j] == κ
                    @constraint(model, g[m,i,κ] ≤ 1 - x[i])
                end
            end
        end
        for j=B[m]
            if r[j] == κ
                @constraint(model, g[m,i,κ] ≤ x[j])
            end
        end
    end

    # Add objective function corresponding to method
    @expression(model, worst_case_expected_revenue, sum(sum(v[m][i]*α[m,i] for i ∈ S[m]) for m=1:M))
    @expression(model, best_case_expected_revenue, sum(κ*g[m,i,κ] for (m,i,κ) ∈ V if m == M))

    if method == 0
        @objective(model, Max, best_case_expected_revenue)
    elseif method == 1
        @objective(model, Max, worst_case_expected_revenue)
    elseif method == 2
        @objective(model, Max, best_case_expected_revenue)
        @constraint(model, worst_case_expected_revenue  ≥ lower_bound)
    end
    
    # Solve optimization problem
    optimize!(model)

    # Extract the optimal solution
    optimal_assortment = Vector{Int64}()
    for i=0:n
        if value(x[i]) > 0.5
            push!(optimal_assortment,i)
        end
    end
    return objective_value(model), optimal_assortment, value(worst_case_expected_revenue)
end



function main()

    # Create a TimerOutput, which we will use to keep track of 
    # the timing and memory usage of various lines of code
    Random.seed!(3)

    N_bar = 0
    E_bar = 0
    r = 0
    sales = 0
    
    # Create the output file
    f_out = open("../data/nested_R1_pareto.csv", "w")

    # Add the headers to the output file
    row_string = string("iter", ",",                      # Iteration that was run

                        "max_previous_assortment", ",",  # Expected revenue of best
														  # previously-offered assortment

    					"optimal_worst_case_revenue", ",",# optimal objective value of
                                                          # robust optimization problem

                        "fraction", ",",                  # fraction of optimal worst-case revenue
    					"worst_case_revenue", ",",		  # Worst-case expected revenue

                        "best_case_revenue", ",",         # Best-case expected revenue

    					"n", ",",                         # Number of products
                        "rev_ordered"                     # true if past assortments are revenue ordered
						)              

    row_string = string(row_string, "\n");
    print(f_out, row_string);
    flush(f_out);

    for iter=1:10
       
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
            Σ,λ_true = GetRandomDistribution(n,80)
   
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
            obj_val_worst_case,optimal_assortment_worst_case,_ = Pareto_Nested(S,sales,r,1)
            
            for fraction in 0.0:0.01:1.0
                
                # Solve the pareto robust optimization problem
                obj_val_pareto,optimal_assortment_pareto, robust_obj_val = Pareto_Nested(S,sales,r,2,obj_val_worst_case*fraction)

                # Compute the worst-case expected revenue for the new assortment
                worst_case,_,_ = Pareto_Nested(S,sales,r,1,0,[],optimal_assortment_pareto)
                
                # Print to file
                row_string = string(iter, ",",
                                    best_past_assortment_obj_val,",",
                                    obj_val_worst_case,",",
                                    fraction,",",
                                    worst_case,",",
                                    obj_val_pareto,",",
                                    n, ",",
                                    setting)
    
                row_string = string(row_string, "\n");
                print(f_out, row_string);
                flush(f_out);
            end
        end
    end
end      
main()
println("") 
