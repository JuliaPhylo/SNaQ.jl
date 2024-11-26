using Graphs, PhyloNetworks

# snaq!(tre0, df, restriction_set(max_level=3, require_galled_tree=true))

function restriction_set(; max_level::Real=Inf, require_galled_tree::Bool=false, require_galled_network::Bool=false)

    restrictions = Vector{Function}()
    if max_level < Inf
        push!(restrictions, restrict_maximum_level(max_level))
    end

    if require_galled_tree
        push!(restrictions, restrict_galled_tree())
    elseif require_galled_network
        push!(restrictions, restrict_galled_network())
    end
    return (net) -> all(F(net) for F in restrictions)

end

restrict_maximum_level(level::Int) = (net) -> get_network_level(net) <= level
restrict_galled_tree() = (net) -> is_galled_tree(net)
restrict_galled_network() = (net) -> is_galled_network(net)


function meets_constraints(net::HybridNetwork, constraints::Vector{Function})
    for F in constraints
        if !F(net) return false end
    end
    return true
end


"""
Converts the tree/network `net` into a SimpleGraph to leverage already
implemented biconnected component algorithms in Graphs.jl.
"""
function Graph(net::HybridNetwork)
    graph = SimpleGraph(net.numnodes)
    nodemap = Dict{Node, Int64}(node => idx for (idx, node) in enumerate(net.node))
    for edge in net.edge
        enode1 = edge.node[1]
        enode2 = edge.node[2]
        if haskey(nodemap, enode1) && haskey(nodemap, enode2)
            add_edge!(graph, nodemap[enode1], nodemap[enode2])
        end
    end

    return graph
end



function get_network_level(net::HybridNetwork)
   
    G = Graph(net)
    bi_comp = biconnected_components(G)
    max_level = 0
    
    for component in bi_comp
        component = unique(reduce(vcat, [[edge.dst, edge.src] for edge in component]))
        iter_level = 0
        for node in component
            if net.node[node].hybrid iter_level += 1 end
        end
        max_level = max(max_level, iter_level)
    end

    return max_level

end


function is_galled_tree(net::HybridNetwork)

    G = Graph(net)
    bi_comp = [comp for comp in biconnected_components(G) if length(comp) > 1]

    v_counts = Dict{Int, Int}()
    for component in bi_comp
        component = reduce(vcat, [[edge.dst, edge.src] for edge in component])
        for node in component
            if !haskey(v_counts, node)
                v_counts[node] = 1
            else
                v_counts[node] += 1
                if v_counts[node] > 2 return false end
            end
        end
    end
    return true

end


function is_galled_network(net::HybridNetwork)

    G = Graph(net)
    bi_comp = [comp for comp in biconnected_components(G) if length(comp) > 1]

    v_counts = Dict{Int, Int}()
    for component in bi_comp
        component = reduce(vcat, [[edge.dst, edge.src] for edge in component])
        for node in component
            if net.node[node].hybrid
                if !haskey(v_counts, node)
                    v_counts[node] = 1
                else
                    v_counts[node] += 1
                    if v_counts[node] > 2 return false end
                end
            end
        end
    end
    return true

end

