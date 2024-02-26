using Graphs
using MetaGraphs
using TimerOutputs

# The following function takes a tuple and a collection of past assortments,
# and outputs a directed acyclic graph in which there is an edge from
# product i to product j if j != i and there exists an assortment m
# such that tuple[m] = i and j in S[m]
# The graph also outputs, for each vertex, the list of vertices
# with an path to that vertex
function ConstructG(assortments, products, tuple)

    # Convert the products and the assortments in products into 
    # strings (because that is needed for doing lookups by
    # product in the graph library)
    products_are_ints = typeof(products[1])  == Int64
    if products_are_ints 
        assortments = [[string(product) for product in assortment] for assortment in assortments]
        products = [string(product) for product in products]
    end

    # Create the graph
    G = MetaDiGraph(length(products))

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


    # Add edges
    for (m,assortment) in enumerate(assortments)
        for i in assortment
            if i != tuple[m]
                MetaGraphs.add_edge!(G, product_ids[i], product_ids[string(tuple[m])])
            end
        end
    end

	# For each vertex i, find the set of vertices that can reach
	# the vertex i_m
    can_reach = Dict()
    for m=1:length(tuple)
        i_m = string(tuple[m])
        can_reach[parse(Int64,i_m)] = Set()
        for i in products
            if i != i_m && has_path(G, product_ids[i], product_ids[i_m])
                push!(can_reach[parse(Int64,i_m)], parse(Int64,i))
            end
        end
    end

	return G, can_reach, product_ids
end



function _GetTuples(G, assortments, product_ids)

	# Base case
	if length(assortments) == 0
        return Set{Tuple}()
	end

	# Create new set of assortments
	tuples = Set{Tuple}()

	# Recursive case
	current_assortment = assortments[1]
	remaining_assortments = assortments[2:length(assortments)]

	# Iterate over the products in current_assortment
	for i in current_assortment
		
		# Find edges to add to the graph
		# Each edge is of the form (j,i), so we only need to
		# keep track of j
		edges_to_add = Set{Any}()

		# Append edges
        for j in current_assortment
            if j != i && !MetaGraphs.has_edge(G, product_ids[j], product_ids[i])
				push!(edges_to_add, j)
                MetaGraphs.add_edge!(G, product_ids[j], product_ids[i])
            end
        end

		# Continue if acyclic
		if !is_cyclic(G)
			if length(remaining_assortments) == 0
				push!(tuples, (i,))
			else
    			sub_tuples = _GetTuples(G, remaining_assortments, product_ids)
    			for sub_tuple in sub_tuples
    				push!(tuples, (i, sub_tuple...))
    			end
			end
		end

		# Remove edges
		for j in edges_to_add
			MetaGraphs.rem_edge!(G, product_ids[j], product_ids[i])
		end
	end

	return tuples
end
function GenerateL(assortments, products)
    # Create the graph
    G = MetaDiGraph(length(products))

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

	# Return the tuples
	return _GetTuples(G, assortments, product_ids)
end

