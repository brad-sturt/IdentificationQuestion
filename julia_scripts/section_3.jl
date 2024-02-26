using Combinatorics
using JuMP
using Gurobi
using Random

include("utils.jl")

function main()
 
    # Get the gurobi environment
    gurobi_env = Gurobi.Env()
 

    ######################################
    # Construct the problem instance
    ######################################
       
    # Specify the number of products
    n = 4

    # Specify the past assortments
    past_assortments = [[0,2,3,4], [0,1,2,4]]

    # Specify the revenues
    revenues = [10, 20, 30, 100]
    r = Dict(i => revenues[i] for i=1:n)
    r[0] = 0

    # Specify the historical sales data
    v = Dict()
    v[2] = Dict(0 => 0.3, 1 => 0.3, 2 => 0.1, 4 => 0.3)
    v[1] = Dict(0 => 0.3, 2 => 0.3, 3 => 0.3, 4 => 0.1)


    ######################################
    # Do the setup work
    ######################################

    # Create the auxiliary sets of products
    S1_minus_S2 = setdiff(Set(past_assortments[1]), Set(past_assortments[2]))
    S2_minus_S1 = setdiff(Set(past_assortments[2]), Set(past_assortments[1]))
    S1_intersect_S2 = intersect(past_assortments[1],past_assortments[2])


    # Construct the set of all possible rankings
    Σ = []
    for perm in permutations(0:n)
        push!(Σ, Dict(i-1 => perm[i] for i=1:(n+1)))
    end

    # Compute the expected revenues for each of the past assortments
    revenues_past_assortments = []
    for m=1:2
        rev = 0
        for i in past_assortments[m]
            rev += v[m][i]*r[i]
        end
        push!(revenues_past_assortments, rev)
    end

    # Compute the expected revenue for the best past assortment
    max_prev_assortment_revenue,max_prev_assortment =  findmax(revenues_past_assortments)

    # Display the expected revenue of the past assortments
    println("--------------------------------------------")
    println("Expected revenue for past assortment 1: ", revenues_past_assortments[1])
    println("Expected revenue for past assortment 2: ", revenues_past_assortments[2])
    println("Expected revenue for best past assortment: ", max_prev_assortment_revenue)

    # Sleep for 3 seconds to see the output
    sleep(1)


    ######################################
    # Compute the worst-case expected revenue
    # for all possible assortments that
    # contain products 0 and 4
    #
    # These numbers are presented in Appendix E.4
    ######################################

    println("--------------------------------------------")

    # Iterate over all of the "possible" assortments, denoted by S_bar
    for S_bar in combinations(0:n)
        
        # If the assortment does not contain 0 and 4, continue
        # to the next assortment
        if 0 ∉ S_bar || 4 ∉ S_bar
            continue
        end

        # Compute the ρ[i1,i2] values for this assortment
        if length(setdiff(intersect(S_bar,past_assortments[2]),past_assortments[1])) != 0
            min_S_intersect_S2_minus_S1 = minimum([r[j] for j in setdiff(intersect(S_bar,past_assortments[2]),past_assortments[1])])
        else
            min_S_intersect_S2_minus_S1 = Inf
        end
        if length(setdiff(intersect(S_bar,past_assortments[1]),past_assortments[2])) != 0
            min_S_intersect_S1_minus_S2 = minimum([r[j] for j in setdiff(intersect(S_bar,past_assortments[1]),past_assortments[2])])
        else
            min_S_intersect_S1_minus_S2 = Inf 
        end
        ρ = Dict()
        for i in S1_intersect_S2
            if i in S_bar
                ρ[i,i] = r[i]
            else
                ρ[i,i] = 0
            end
        end
        for i1 in S1_intersect_S2, i2 in S2_minus_S1
            if i2 in S_bar
                ρ[i1,i2] = r[i2]
            elseif i1 in S_bar
                ρ[i1,i2] = min(r[i1], min_S_intersect_S2_minus_S1)
            else
                ρ[i1,i2] = 0
            end
        end
        for i1 in S1_minus_S2, i2 in S1_intersect_S2
            if i1 in S_bar
                ρ[i1,i2] = r[i1]
            elseif i2 in S_bar
                ρ[i1,i2] = min(r[i2], min_S_intersect_S1_minus_S2)
            else
                ρ[i1,i2] = 0
            end
        end
        for i1 in S1_minus_S2, i2 in S2_minus_S1
            if i1 in S_bar && i2 in S_bar
                ρ[i1,i2] = min(r[i1],r[i2])
            elseif i1 in S_bar && i2 ∉ S_bar
                ρ[i1,i2] = min(r[i1], min_S_intersect_S2_minus_S1)
            elseif i1 ∉ S_bar && i2 in S_bar
                ρ[i1,i2] = min(r[i2], min_S_intersect_S1_minus_S2)
            else
                ρ[i1,i2] = 0
            end
        end

        # Construct the robust optimization problem
        model = direct_model(Gurobi.Optimizer(gurobi_env))::JuMP.Model
        set_silent(model)
        @variable(model, λ[i1=past_assortments[1],i2=past_assortments[2]] ≥ 0)

        for i1 in past_assortments[1], i2 in past_assortments[2]
            if !(i1 in S1_minus_S2 && i2 in S2_minus_S1) && !(i1 in S1_minus_S2 && i2 in S1_intersect_S2) && !(i1 in S1_intersect_S2 && i2 in S2_minus_S1)
                @constraint(model, λ[i1,i2] == 0)
            end
        end
        @constraint(model, Block1[i1 in S1_minus_S2], sum(λ[i1,i2] for i2 in past_assortments[2]) == v[1][i1])
        @constraint(model, Block2[i2 in S2_minus_S1], sum(λ[i1,i2] for i1 in past_assortments[1]) == v[2][i2])
        @constraint(model, Block3[i in S1_intersect_S2], sum(λ[i1,i] for i1 in S1_minus_S2) - sum(λ[i,i2] for i2 in S2_minus_S1) == v[2][i] - v[1][i])
        @objective(model, Min, sum(sum(ρ[i1,i2]*λ[i1,i2] for i1 in S1_minus_S2) for i2 in S1_intersect_S2) + 
                               sum(sum(ρ[i1,i2]*λ[i1,i2] for i1 in S1_minus_S2) for i2 in S2_minus_S1) + 
                               sum(sum((ρ[i1,i2] - ρ[i1,i1])*λ[i1,i2] for i1 in S1_intersect_S2) for i2 in S2_minus_S1))

        # Solve the model
        optimize!(model)

        # Get the optimal objective value
        robust_revenue = objective_value(model) + sum(ρ[i,i]*v[1][i] for i in S1_intersect_S2)

        # Print the assortment and corresponding worst-case expected revenue
        println(sort(collect(S_bar)), "\t", robust_revenue)
    end


    #########################################################
    # Repeatedly perform estimate-then-optimize to
    # generate many estimates of the ranking-based choice
    # model (and thus obtain various different new assortments
    # to investigate) until we obtain the assortment {0,4}
    #
    # The specific ranking-based choice models from Appendix 
    # E.2 and E.3 were found by the output of the following code
    #########################################################

    println("--------------------------------------------")

    # Loop until estimate-then-optimize outputs the assormtent {0,4}
    while true
    
        ###################################################
     	# Find assortment using estimate-then-optimize
        ###################################################

        # Compute an estimate of the ranking-based choice model
        λ_est = EstimateDistribution(past_assortments,v,Σ,gurobi_env)

        # Compute the assortment which maximizes the expected revenue under the
        # estimated ranking-based choice model
        assortment_est, revenue_est = GetOptimalAssortment(λ_est,n,r,Σ,gurobi_env)
        assortment_est = sort(collect(assortment_est))

        # If we did not obtain the assortment {0,4}, go back to beginning of loop
        if assortment_est != [0,4]
            continue
        end
        
        ###################################################
        # Show the parameters of the ranking-based choice
        # model that generated this assortment.
        ##################################################
        
        # Obtain the rankings which have nonzero probability under
        # the estimated ranking-based choice model
        estimated_dist_indices = [i for i=1:factorial(n+1) if λ_est[i] > 1e-8]

        # Display the parameters and rankings
        for i=estimated_dist_indices
            @show Σ[i], λ_est[i]
        end


        ###################################################
     	# Show the parameters for a ranking-based
        # choice model that generates the worst-case expected 
        # revenue of the assortment obtained using estimate-
        # then-optimize
        ###################################################
       
        println("--------------------------------------------")

        # Get the worst-case ranking-based choice model
        revenue_wc, λ_wc = EvaluateAssortment(assortment_est,r,past_assortments,v,n,Σ,gurobi_env,false)

        # Obtain the rankings which have nonzero probability under
        # the worst-case ranking-based choice model
        wc_dist_indices = [i for i=1:factorial(n+1) if λ_wc[i] > 1e-8]

        # Display the parameters and rankings
        for i=wc_dist_indices
            @show Σ[i], λ_wc[i]
        end

        break
    end
end
