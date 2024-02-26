
###########################################################
# Construct a dictionary, denoted by `ρ`, which maps
# each tuple to the worst-case revenue of a customer
# whose preferences correspond to the tuple
###########################################################

function Construct_ρ(assortment, L, tuple_to_can_reach,  r)
    ρ = Dict()
    for tuple in L
        ρ[tuple] = Inf
        can_reach = tuple_to_can_reach[tuple]
        if tuple == (0,0)
        end
        for i in assortment
            # Create flag to check whether it can be added
            can_be_added = true
            for i_m in tuple
                if i_m in assortment
                    if i in can_reach[i_m]
                        can_be_added = false
                        break
                    end
                end
            end
            if can_be_added
                ρ[tuple] = min(ρ[tuple], r[i])
            end
        end
    end
    return ρ
end

###########################################################
# Construct a dictionary, denoted by `ξ`, which maps
# each tuple to the best-case revenue of a customer
# whose preferences correspond to the tuple
###########################################################

function Construct_ξ(assortment, L, tuple_to_can_reach,  r)
    ξ = Dict()
    for tuple in L
        ξ[tuple] = -Inf
        can_reach = tuple_to_can_reach[tuple]
        for i in assortment
            # Create flag to check whether it can be added
            can_be_added = true
            for i_m in tuple
                if i_m in assortment
                    if i in can_reach[i_m]
                        can_be_added = false
                        break
                    end
                end
            end
            if can_be_added
                ξ[tuple] = max(ξ[tuple], r[i])
            end
        end
    end
    return ξ
end
