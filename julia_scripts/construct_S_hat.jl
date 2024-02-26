using Graphs
using MetaGraphs
using TimerOutputs

function ConstructG_S_hat(assortments, products, r)
    
    # Create the graph
    G = MetaDiGraph(length(products))

    # Convert the products and the assortments in products into 
    # strings (because that is needed for doing lookups by
    # product in the graph library)
    #products_are_ints = typeof(products[1]) == Int64
    #if products_are_ints 
        assortments = [[string(product) for product in assortment] for assortment in assortments]
        products = [string(product) for product in products]
        r= Dict(string(product) => r[product] for product in keys(r))
    #end

    # Create a dictionary that maps each product to an initial product
    # id starting at 1
    product_ids = Dict{Any,Int64}()
    for (product_id,product) in enumerate(products)
        product_ids[product] = product_id
    end

    # Add the product names as attributes to each vertex in the digraph
    for product in products
        set_prop!(G,product_ids[product],:name,product)
    end

    # Set the indexing for G by name
    set_indexing_prop!(G, :name)

    # Construct the set M_i of previously-offered assortments which offered 
    # product i
    M = Dict{Any,Set{Any}}()
    for i in products
        M[i] = Set{Any}()
    end
    for assortment in assortments
        for i in assortment
            push!(M[i], assortment)
        end 
    end

    # For each pair of products, add an edge if necessary
    for i_star in products
        for i in products
            if r[i_star] < r[i] && intersect(M[i_star], M[i]) == M[i_star]
                MetaGraphs.add_edge!(G, product_ids[i_star], product_ids[i])
            end
        end
    end

    return G
end

# Note that the following function will break if the product names are not
# strings. 
function RecursiveStep(G)

    # Base case
    if MetaGraphs.nv(G) == 0
        A = Set{Set{Any}}()
        push!(A, Set{Any}())
        return A
    end
   
    # Get the first vertex in G
    i = G[1,:name]

    # Get the list of incoming vertices
    incoming_to_i = [G[ℓ,:name] for ℓ in inneighbors(G,G[i,:name])]

    # Get the list of outgoing vertices
    outgoing_from_i = [G[ℓ,:name] for ℓ in outneighbors(G,G[i,:name])]

    # Create a copy of G in which the vertices {i} ∪ {ℓ: (ℓ,i) ∈ edges(G)} 
    # and the incoming and outgoing edges of those vertices are removed.
    G_prime = deepcopy(G)

    # Remove i from G_prime
    rem_vertex!(G_prime,G_prime[i,:name])

    # Remove incoming vertices from G_prime
    for ℓ in incoming_to_i
        rem_vertex!(G_prime,G_prime[ℓ,:name])
    end

    # Do recursion on G_prime
    A_prime = RecursiveStep(G_prime)
 
    # Create a copy of G in which the vertices {i} ∪ {ℓ: (i,ℓ) ∈ edges(G)} 
    # and the incoming and outgoing edges of those vertices are removed.
    G_prime_prime = deepcopy(G)

    # Remove i from G_prime_prime
    rem_vertex!(G_prime_prime,G_prime_prime[i,:name])

    # Remove incoming vertices from G_prime_prime
    for ℓ in outgoing_from_i
        rem_vertex!(G_prime_prime,G_prime_prime[ℓ,:name])
    end

    # Do recursion on G_prime
    A_prime_prime = RecursiveStep(G_prime_prime)

    # Create the collection A_prime_prime_prime
    A_prime_prime_prime = Set{Set{Any}}()
    for S in A_prime_prime
        S_new = Set{Any}()
        push!(S_new, i)
        for ℓ in outgoing_from_i
            push!(S_new, ℓ)
        end
        S_new = union(S_new, S)
        push!(A_prime_prime_prime, S_new)
    end
  
    # Return
    return union(A_prime, A_prime_prime_prime)
end

function ConstructS_hat(assortments, products, r)

    assortments_simplified = assortments
    r_simplified = r
    products_simplified = products
    
    # Construct the directed acyclic graph
    G = ConstructG_S_hat(assortments_simplified, products_simplified, r_simplified)

    # Compute the collection of assortments
    S_hat_temp = RecursiveStep(G)

    # Add the removed products to each of the assortments in S_hat
    S_hat = Set{Set{Any}}()
    for S_temp in S_hat_temp
        S = Set([parse(Int64,product) for product in S_temp])
        if 0 ∉ S #|| length(products)-1 ∉ S
            continue
        end
        #=
        for product in products_that_appear_in_every_assortment
            push!(S, product)
        end
        =#
        push!(S_hat, S)
    end

    # Return 
    return S_hat
    
end

function ConstructG_S_hat_optimistic(assortments, products, r)
    
    # Create the graph
    G = MetaDiGraph(length(products)-1)

    # Convert the products and the assortments in products into 
    # strings (because that is needed for doing lookups by
    # product in the graph library)
    assortments = [[string(product) for product in assortment if product != 0] for assortment in assortments]
    products = [string(product) for product in products if product != 0]
    r= Dict(string(product) => r[product] for product in keys(r))

    # Create a dictionary that maps each product to an initial product
    # id starting at 1
    product_ids = Dict{Any,Int64}()
    for (product_id,product) in enumerate(products)
        product_ids[product] = product_id
    end

    # Add the product names as attributes to each vertex in the digraph
    for product in products
        if product != "0"
            set_prop!(G,product_ids[product],:name,product)
        end
    end

    # Set the indexing for G by name
    set_indexing_prop!(G, :name)

    # Construct the set M_i of previously-offered assortments which offered 
    # product i
    M = Dict{Any,Set{Any}}()
    for i in products
        M[i] = Set{Any}()
    end
    for assortment in assortments
        for i in assortment
            push!(M[i], assortment)
        end 
    end

    # For each pair of products, add an edge if necessary
    for i_star in products
        for i in products
            if i_star == "0" || i == "0"
                continue
            end
            if r[i_star] < r[i] && intersect(M[i_star], M[i]) == M[i]
                MetaGraphs.add_edge!(G, product_ids[i_star], product_ids[i])
            end
        end
    end

    return G
end

# Note that the following function will break if the product names are not
# strings. 
function RecursiveStep_optimistic(G)

    # Base case
    if MetaGraphs.nv(G) == 0
        A = Set{Set{Any}}()
        push!(A, Set{Any}())
        return A
    end
   
    # Get the first vertex in G
    i = G[1,:name]

    # Get the list of incoming vertices
    incoming_to_i = [G[ℓ,:name] for ℓ in inneighbors(G,G[i,:name])]

    # Get the list of outgoing vertices
    outgoing_from_i = [G[ℓ,:name] for ℓ in outneighbors(G,G[i,:name])]

    # Create a copy of G in which the vertex {i}  
    # and the incoming and outgoing edge are removed.
    G_prime = deepcopy(G)

    # Remove i from G_prime
    rem_vertex!(G_prime,G_prime[i,:name])

    # Remove incoming vertices from G_prime

    # Do recursion on G_prime
    A_prime = RecursiveStep_optimistic(G_prime)
 
    # Create a copy of G in which the vertices {i} ∪ {ℓ: (i,ℓ) ∈ edges(G)} 
    # and the incoming and outgoing edges of those vertices are removed.
    G_prime_prime = deepcopy(G)

    # Remove i from G_prime_prime
    rem_vertex!(G_prime_prime,G_prime_prime[i,:name])

    # Remove incoming vertices from G_prime_prime
    for ℓ in outgoing_from_i
        rem_vertex!(G_prime_prime,G_prime_prime[ℓ,:name])
    end
    for ℓ in incoming_to_i
        rem_vertex!(G_prime_prime,G_prime_prime[ℓ,:name])
    end

    # Do recursion on G_prime
    A_prime_prime = RecursiveStep_optimistic(G_prime_prime)

    # Create the collection A_prime_prime_prime
    A_prime_prime_prime = Set{Set{Any}}()
    for S in A_prime_prime
        S_new = Set{Any}()
        push!(S_new, i)
        S_new = union(S_new, S)
        push!(A_prime_prime_prime, S_new)
    end
  
    # Return
    return union(A_prime, A_prime_prime_prime)
end

function ConstructS_hat_optimistic(assortments, products, r)
   
    # Construct the directed acyclic graph
    G = ConstructG_S_hat_optimistic(assortments, products, r)

    # Compute the collection of assortments
    S_hat_temp = RecursiveStep_optimistic(G)

    # Add the removed products to each of the assortments in S_hat
    S_hat = Set{Set{Any}}()
    for S_temp in S_hat_temp
        S = Set([parse(Int64,product) for product in S_temp])
        #if 0 ∉ S #|| length(products)-1 ∉ S
        if length(products)-1 ∉ S
            continue
        end
        push!(S,0)
        #=
        for product in products_that_appear_in_every_assortment
            push!(S, product)
        end
        =#
        push!(S_hat, S)
    end

    # Return 

    return S_hat
    
end

