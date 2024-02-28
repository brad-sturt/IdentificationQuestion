using Combinatorics
using JuMP
using Gurobi
using Random

include("utils.jl")

function main()
 
    ###################################################
 	# Initialize the output file, do other 
    # miscellaneous initializations
	###################################################
   
    # Set seed
    Random.seed!(123)
    
    # Get the gurobi environment
    gurobi_env = Gurobi.Env()
    
    # Create the output file
    f_out = open("../data/appendix_A2_2.csv", "w")
    
    
    # Add the headers to the output file
    row_string = string("iter", ",",  # Iteration that was run

                        "n", ",",     # Number of products

                        "K", ",",     # Number of rankings with nonzero support
                        
                        "speed"       # Computation time in seconds
                        )
    row_string = string(row_string, "\n");
    print(f_out, row_string);
    flush(f_out);

    
    ###################################################
 	# Perform each iteration for each value of n
	###################################################

    # Set the number of iterations, choice of K, and choices of
    # the number of products n
    num_iterations = 100
    Ks = [1000,]
    ns = 10:2:100

    for iter=1:num_iterations, K in Ks, n in ns


        ###################################################
     	# Initialize the problem instance
    	###################################################
     
        # Generate random revenues for the products
        revenues = rand(n)
        revenues = sort(revenues)
        r = Dict(i => revenues[i] for i=1:n)
        r[0] = 0
    
        # Generate two random past assortments
        indices = rand(1:3,n)  # 1 -> past_assortments[1]
                               # 2 -> past_assortments[2]
                               # 3 -> past_assortments[1] and past_assortments[2]
        past_assortments = [[0],[0]]
        S1_minus_S2 = Set{Int64}()
        S2_minus_S1 = Set{Int64}()
        S1_intersect_S2 = Set([0,n])
        for i=1:(n-1)
            if indices[i] == 1
                push!(past_assortments[1],i)
                push!(S1_minus_S2,i)
            elseif indices[i] == 2
                push!(past_assortments[2],i)
                push!(S2_minus_S1,i)
            else
                push!(past_assortments[1],i)
                push!(past_assortments[2],i)
                push!(S1_intersect_S2,i)
            end
        end
        push!(past_assortments[1],n)
        push!(past_assortments[2],n)


        # Generate a random probability distribution
        Σ, λ_base = GetRandomDistribution(n,K)

        # Compute the historical sales data
        v = GetFrequencies(λ_base,past_assortments,Σ)

        time_elapsed =  @elapsed begin

            # Get the revenues for each of the previously-offered assortments
            revenues_past_assortments = []
            for m=1:2
                rev = 0
                for i in past_assortments[m]
                    rev += v[m][i]*r[i]
                end
                push!(revenues_past_assortments, rev)
            end
            max_prev_assortment_revenue,max_prev_assortment =  findmax(revenues_past_assortments)
             
            
    
            ###################################
            # Get optimal robust assortment
            ###################################
    
            # Create the temporary variables
            best_robust_assortment_revenue = 0
            best_robust_assortment = []
    
    
            # Iterate over all of the "possible" assortments
            for i1_temp=union(S1_minus_S2,Set(Inf)),i2_temp=union(S2_minus_S1,Set(Inf))
    
                # Construct the assortment
                S_bar = Set()
                for j in S1_minus_S2
                    if j >= i1_temp
                        push!(S_bar, j)
                    end
                end
                for j in S2_minus_S1
                    if j >= i2_temp
                        push!(S_bar,j)
                    end
                end
                for j in S1_intersect_S2
                    push!(S_bar,j)
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

                # Following line added in second revision
                @constraint(model, Block4[i in S1_intersect_S2], sum(λ[i,i2] for i2 in S2_minus_S1) ≤ v[1][i])

                @objective(model, Min, sum(sum(ρ[i1,i2]*λ[i1,i2] for i1 in S1_minus_S2) for i2 in S1_intersect_S2) + 
                                       sum(sum(ρ[i1,i2]*λ[i1,i2] for i1 in S1_minus_S2) for i2 in S2_minus_S1) + 
                                       sum(sum((ρ[i1,i2] - ρ[i1,i1])*λ[i1,i2] for i1 in S1_intersect_S2) for i2 in S2_minus_S1))
    
                # Solve the model
                optimize!(model)
    
                # Get the optimal objective value
                robust_revenue = objective_value(model) + sum(ρ[i,i]*v[1][i] for i in S1_intersect_S2)
    
    
    
    
                # Continue if necessary
                if robust_revenue ≤ best_robust_assortment_revenue
                    continue
                end
    
                # Update if necessary
                best_robust_assortment_revenue = robust_revenue
                best_robust_assortment = sort(collect(S_bar))
    
            end
            #_,revenue_wc = Robust(assortment_est,r,S,v,A,n,Σ_all,gurobi_env)
            #println(best_robust_assortment_revenue, ", ", max_prev_assortment_revenue)
        end    
        @show time_elapsed, n
    
        row_string = string(iter, ",",
                            n, ",",
                            K, ",",
                            time_elapsed
                            )
        row_string = string(row_string, "\n");
        print(f_out, row_string);
        flush(f_out);
    end
end
