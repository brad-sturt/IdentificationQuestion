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



function SolveBestCase_Nested(S,v,r)

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
        if r[i] > r[i+1]
            error("The revenues for the products are not sorted")
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
    @expression(model, best_case_expected_revenue, sum(κ*g[m,i,κ] for (m,i,κ) ∈ V if m == M))

    @objective(model, Max, best_case_expected_revenue)
    
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


function SolveWorstCase_Nested(S,v,r)

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
        if r[i] > r[i+1]
            error("The revenues for the products are not sorted")
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

    @expression(model, worst_case_expected_revenue, sum(sum(v[m][i]*α[m,i] for i ∈ S[m]) for m=1:M))

    @objective(model, Max, worst_case_expected_revenue)

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
        if r[i] > r[i+1]
            error("The revenues for the products are not sorted")
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

