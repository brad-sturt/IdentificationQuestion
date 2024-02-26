using JuMP
using Gurobi
using Combinatorics

using Printf
include("construct_rho.jl")
include("construct_S_hat.jl")

function SolveWorstCase(new_assortment, assortments, products, ρ, v,L, gurobi_env)
    M = length(assortments)
    
    # Create the model
    model = direct_model(Gurobi.Optimizer(gurobi_env))::JuMP.Model
    set_silent(model)

    # Variables
    @variable(model, λ[L] ≥ 0)
    
    # Constraints
    # Creat mapping from each (m,i) to relevant tuples
    tuple_mapping = Dict()
    for m=1:M, i=assortments[m]
        tuple_mapping[m,i] = Set()
    end
    for tuple in L
        for m=1:M
            push!(tuple_mapping[m,tuple[m]], tuple)
        end
    end
    @constraint(model, Equality[m=1:M,i=assortments[m]], sum(λ[tuple] for tuple in tuple_mapping[m,i]) == v[m][i])
    
    # Objective
    @objective(model, Min, sum(ρ[tuple]*λ[tuple] for tuple in L))

    # Solve the model
    optimize!(model)
    return objective_value(model)
end

function SolveRobustOptimizationProblem(assortments, products, v, r, gurobi_env)

 	# Construct L
	L = GenerateL(assortments,products)

	# Construct directed cyclic graph for each tuple
	tuple_to_dag = Dict()
	tuple_to_product_ids = Dict()
	tuple_to_can_reach = Dict()
	for tuple in L
		tuple_to_dag[tuple],tuple_to_can_reach[tuple], tuple_to_product_ids[tuple] = ConstructG(assortments,products, tuple)
	end

    # Construct S_hat
    S_hat = ConstructS_hat(assortments, products, r)
    S_hat_cleaned = sort([sort(collect(S)) for S in S_hat if sort(collect(S)) ∉ assortments])

    # Keep track of best objective value 
    best_obj_val = -Inf
    best_assortment = []
    assortment_to_obj_val = Dict()

    # Iterate over the assortments in S_hat
    for assortment in S_hat

        # Convert assortment to string
        ρ = Construct_ρ(assortment, L, tuple_to_can_reach, r)
    
        obj_val =  SolveWorstCase(assortment, assortments, products, ρ, v, L, gurobi_env)
        assortment_to_obj_val[sort(collect(assortment))] = obj_val
        if obj_val > best_obj_val
            best_obj_val = obj_val
            best_assortment = assortment
        end
    end
    S_hat_clean = sort([assortment for assortment in keys(assortment_to_obj_val)])

    return best_obj_val, best_assortment, S_hat
end

function EvaluateWorstAndBestCases(new_assortments, assortments, products, v, r, gurobi_env)
 	# Construct L
	L = GenerateL(assortments,products)

	# Construct directed cyclic graph for each tuple
	tuple_to_dag = Dict()
	tuple_to_product_ids = Dict()
	tuple_to_can_reach = Dict()
	for tuple in L
		tuple_to_dag[tuple],tuple_to_can_reach[tuple], tuple_to_product_ids[tuple] = ConstructG(assortments,products, tuple)
	end

    # Keep track of best objective value 
    assortment_to_ro_obj_val = Dict()
    assortment_to_oo_obj_val = Dict()

    # Iterate over the assortments in S_hat
    for assortment in new_assortments

        # Compute objective function for worst-case
        ρ = Construct_ρ(assortment, L, tuple_to_can_reach, r)
 
        # Compute objective function for best-case
        ξ = Construct_ξ(assortment, L, tuple_to_can_reach, r)
 
        # Compute the worst-case and best-case objective value
        assortment_to_ro_obj_val[sort(collect(assortment))] = SolveWorstCase(assortment, assortments, products, ρ, v, L, gurobi_env)
        assortment_to_oo_obj_val[sort(collect(assortment))] = SolveBestCase(assortment, assortments, products, ξ, v, L, gurobi_env)

    end

    return assortment_to_ro_obj_val, assortment_to_oo_obj_val 

end


function SolveBestCase(new_assortment, assortments, products, ξ, v,L, gurobi_env)


    M = length(assortments)
    
    # Create the model
    model = direct_model(Gurobi.Optimizer(gurobi_env))::JuMP.Model
    set_silent(model)

    # Variables
    @variable(model, λ[L] ≥ 0)
    
    # Constraints
    @constraint(model, Equality[m=1:M,i=assortments[m]], sum(λ[tuple] for tuple in L if tuple[m] == i) == v[m][i])

    # Objective
    @objective(model, Max, sum(ξ[tuple]*λ[tuple] for tuple in L))

    # Solve the model
    optimize!(model)
    return objective_value(model)
end


function SolveOptimisticOptimizationProblem(assortments, products, v, r, gurobi_env)

 	
    # Construct L
	L = GenerateL(assortments,products)

	# Construct directed cyclic graph for each tuple
	tuple_to_dag = Dict()
	tuple_to_product_ids = Dict()
	tuple_to_can_reach = Dict()
	for tuple in L
		tuple_to_dag[tuple],tuple_to_can_reach[tuple], tuple_to_product_ids[tuple] = ConstructG(assortments,products, tuple)
	end

    # Construct S_hat
    S_hat_optimistic = ConstructS_hat_optimistic(assortments, products, r)

    # Keep track of best objective value 
    best_obj_val = -Inf
    best_assortment = []

    assortment_to_obj_val = Dict()

    for assortment in S_hat_optimistic
        if 0 ∉ assortment || length(products)-1 ∉ assortment
            #error()
            @show assortment
            continue
        end

        # Convert assortment to string
        ξ = Construct_ξ(assortment, L, tuple_to_can_reach, r)
  
        obj_val = SolveBestCase(assortment, assortments, products, ξ, v, L,gurobi_env)
        #=
        print(assortment, "\n")
        for tuple in L
            print(ρ[tuple], "\t")
        end
        print("\n")
        #@show obj_val, assortment
        =#
        assortment_to_obj_val[sort(collect(assortment))] = obj_val
        if obj_val > best_obj_val
            best_obj_val = obj_val
            best_assortment = assortment
        end
    end
    S_hat_optimistic_clean = sort([assortment for assortment in keys(assortment_to_obj_val)])
    for assortment in S_hat_optimistic_clean
       #println(@sprintf("%-30s %4.5f",assortment, assortment_to_obj_val[assortment]))
        #println(assortment, ",\t", assortment_to_obj_val[assortment])
    end

    return best_obj_val, best_assortment, S_hat_optimistic
end
