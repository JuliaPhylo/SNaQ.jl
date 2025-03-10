include("CF_struct.jl")
include("CF_blocks.jl")
include("CF_recursive_blocks.jl")
using Graphs


function get_reticulate_4taxa_quartet_equations(net::HybridNetwork, taxa::Vector{<:AbstractString}, edge_number_to_idx_map::Dict{Int, Int})

    if net.numhybrids == 0
        quartet_type, int_edges = get_quartet_type_and_internal_edges(net, taxa, edge_number_to_idx_map)
        return RecursiveCFEquation(
            true, int_edges, quartet_type, nothing, Vector{RecursiveCFEquation}([])
        )
    end

    # If >= 1 hybrids, we need to keep recursing
    lowest_H = get_lowest_hybrid(net)
    if n_leaves_below_lowest_hybrid(lowest_H) == 1

        # Remove the minor edge and all of its references in this copy
        div_major = deepcopy(net)
        div_major_H = div_major.hybrid[findfirst(h -> h.name == lowest_H.name && h.number == lowest_H.number, div_major.hybrid)]
        E_minor = getparentedgeminor(div_major_H)
        E_major = getparentedge(div_major_H)
        for node in E_minor.node
            node.edge = [e for e in node.edge if e != E_minor]
        end
        PN.deleteEdge!(div_major, E_minor; part=false)
        PN.removeHybrid!(div_major, getchild(E_major))
        E_major.hybrid = false
        getchild(E_major).hybrid = false

        # Remove the major edge and all of its references in this copy
        div_minor = deepcopy(net)
        div_minor_H = div_minor.hybrid[findfirst(h -> h.name == lowest_H.name && h.number == lowest_H.number, div_minor.hybrid)]
        E_minor = getparentedgeminor(div_minor_H)
        E_major = getparentedge(div_minor_H)
        for node in E_minor.node
            node.edge = [e for e in node.edge if e != E_minor]
        end
        PN.deleteEdge!(div_minor, E_minor; part=false)
        PN.removeHybrid!(div_minor, getchild(E_major))   # only removes its references - does not delete the node
        E_major.hybrid = false
        getchild(E_major).hybrid = false

        return RecursiveCFEquation(
            false, [], 0, lowest_H, [
                get_reticulate_4taxa_quartet_equations(div_minor, taxa, edge_number_to_idx_map),
                get_reticulate_4taxa_quartet_equations(div_major, taxa, edge_number_to_idx_map)
            ]
        )

    else
        error("Implemented 2/4 cases where there are 2 leaves below hybrid so far - need to implement remaining 2 cases.")
        int_edges = get_internal_edges_below_lowest_hybrid(lowest_H)
        leaves_below_H = get_leaves_below_lowest_hybrid(lowest_H)

        # Both take minor:
        # 1. remove internal edges (don't fuse - we need 
        #    edges to still map to the original net)
        # 2. remove major retic edge
        # 3. H and its retics are no longer hybrids
        # 4. each leaf below H is now child to H
        div1 = deepcopy(net)
        div1_int_edges = div1.edge[[findfirst(div1_edge -> div1_edge.number == int_edge.number, div1.edge) for int_edge in int_edges]]
        div1_H = div1.hybrid[findfirst(div1_H -> div1_H.number == lowest_H.number && div1_H.name == lowest_H.name, div1.hybrid)]
        E_minor = getparentedgeminor(div1_H)
        E_major = getparentedge(div1_H)
        for node in E_major.node
            node.edge = [e for e in node.edge if e != E_major]
        end
        PN.deleteEdge!(div1, E_major; part=false)
        PN.removeHybrid!(div1, getchild(E_minor))   # only removes its references - does not delete the node
        E_minor.hybrid = false
        E_minor.ismajor = true
        getchild(E_minor).hybrid = false

        for L in leaves_below_H
            new_leaf = PN.addleaf!(div1, getchild(div1_H), "__$(L.name)", 0.0)
            PN.deleteleaf!(div1, L.name; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            new_leaf.name = L.name
        end

        # Both take major
        # 1. remove internal edges
        # 2. remove minor retic edge
        # 3. H and its retics are no longer hybrids
        # 4. each leaf below H is now child to H
        div2 = deepcopy(net)
        div2_int_edges = div2.edge[[findfirst(div2_edge -> div2_edge.number == int_edge.number, div2.edge) for int_edge in int_edges]]
        div2_H = div2.hybrid[findfirst(div2_H -> div2_H.number == lowest_H.number && div2_H.name == lowest_H.name, div2.hybrid)]
        E_minor = getparentedgeminor(div2_H)
        E_major = getparentedge(div2_H)
        for node in E_minor.node
            node.edge = [e for e in node.edge if e != E_minor]
        end
        PN.deleteEdge!(div2, E_minor; part=false)
        PN.removeHybrid!(div2, getchild(E_major))   # only removes its references - does not delete the node
        E_major.hybrid = false
        getchild(E_major).hybrid = false

        for L in leaves_below_H
            new_leaf = PN.addleaf!(div2, getchild(div2_H), "__$(L.name)", 0.0)
            PN.deleteleaf!(div2, L.name; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            new_leaf.name = L.name
        end

        @error "Only implemented 2/4 parental tree splits for test - still need last 2"
        leaf_names = sort([leaves_below_H[1].name, leaves_below_H[2].name])
        which_quartet = leaf_names[1] == taxa[1] ? (
            leaf_names[2] == taxa[2] ? 1 :
            leaf_names[2] == taxa[3] ? 2 : 3
        ) :
        leaf_names[1] == taxa[2] ? (
            leaf_names[2] == taxa[3] ? 3 : 2
        ) : 1
        return RecursiveCFEquation(
            length(int_edges) > 0, int_edges, which_quartet, lowest_H, [
                get_reticulate_4taxa_quartet_equations(div1, taxa, edge_number_to_idx_map),
                get_reticulate_4taxa_quartet_equations(div2, taxa, edge_number_to_idx_map)
            ]
        )


        # Lower taxa takes minor, higher takes major


    end

end


function get_quartet_type_and_internal_edges(net::HybridNetwork, taxa::Vector{<:AbstractString}, edge_number_to_idx_map::Dict{Int, Int})
    # iterate over the 2^H displayed trees
    G, W = Graph(net; withweights=true, minoredgeweight=Inf)
    for idx in eachindex(W) W[idx] = (W[idx] == Inf) ? Inf : 1.0 end
    node_to_idx = Dict{PN.Node, Int}(node => j for (j, node) in enumerate(net.node))                                            # these two dicts used later for
    edge_to_graph_idxs = Dict{PN.Edge,Tuple{Int,Int}}(e => (node_to_idx[e.node[1]], node_to_idx[e.node[2]]) for e in net.edge)  # easier graph weight adjustment
    
    # Find which edges form the internal edge of the quartet
    taxa1_node_idx = findfirst(n -> n.leaf && n.name == taxa[1], net.node)
    taxa2_node_idx = findfirst(n -> n.leaf && n.name == taxa[2], net.node)
    taxa3_node_idx = findfirst(n -> n.leaf && n.name == taxa[3], net.node)
    taxa4_node_idx = findfirst(n -> n.leaf && n.name == taxa[4], net.node)

    path_12 = a_star(G, taxa1_node_idx, taxa2_node_idx, W)
    path_34 = a_star(G, taxa3_node_idx, taxa4_node_idx, W)
    path_13 = a_star(G, taxa1_node_idx, taxa3_node_idx, W)
    path_24 = a_star(G, taxa2_node_idx, taxa4_node_idx, W)

    # edges here have direction, but we don't want them to, so we set (src,dst) = (minidx,maxidx) so the direction is effectively removed
    for path in [path_12, path_34, path_13, path_24]
        for (j, edge) in enumerate(path)
            path[j] = Graphs.SimpleEdge(min(edge.src, edge.dst), max(edge.src, edge.dst))
        end
    end

    if length(intersect(path_12, path_34)) == 0
        internal_graph_edges = intersect(path_13, path_24)    # 12|34 is displayed, so these paths must cross ONLY on the displayed edge portion
        internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)
        return 1, internal_net_edges
    elseif length(intersect(path_13, path_24)) == 0
        internal_graph_edges = intersect(path_12, path_34)    # 13|24 is displayed, so these paths must cross ONLY on the displayed edge portion
        internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)
        return 2, internal_net_edges
    else
        internal_graph_edges = intersect(path_12, path_34)    # 14|23 is displayed, so these paths must cross ONLY on the displayed edge portion
        internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)
        return 3, internal_net_edges
    end
end


"""
Gets the "lowest" hybrid, i.e. one of potentially multiple hybrids that do not have any other hybrids in their descendants.
Function assumes that extraneous retics have already been removed (i.e. retics on external quartet branches).
"""
function get_lowest_hybrid(net::HybridNetwork)::Node
    if net.numhybrids == 1 return net.hybrid[1] end
    return get_lowest_hybrid_recur(net.hybrid[1])
end

function get_lowest_hybrid_recur(node::Node)
    if node.leaf
        return nothing
    end

    children = getchildren(node)
    for child in children
        child_val = get_lowest_hybrid_recur(child)
        if child_val !== nothing return child_val end
    end

    if node.hybrid
        return node
    else
        return nothing
    end
end


"""
Gets the number of leaves in a quarnet below the lowest hybrid in the quarnet.
Assumes that reticulations on external edges are removed.
"""
function n_leaves_below_lowest_hybrid(H::Node)
    c = getchildren(H)
    while length(c) == 1
        if c[1].leaf
            return 1
        else
            c = getchildren(c[1])
        end
    end
    return 2
end


"""
Assumes that there are 2 leaves below `H` in the quarnet.
"""
function get_internal_edges_below_lowest_hybrid(H::Node)
    internal_edges = Vector{Edge}()
    c = getchildren(H)
    while length(c) == 1
        push!(internal_edges, getparentedge(c[1]))
        c = getchildren(c[1])
    end
    push!(internal_edges, getparentedge(c[1]))
    return internal_edges
end


function get_leaves_below_lowest_hybrid(H::Node)
    queue = Vector{Node}([H])
    leaves = Vector{Node}([])

    while length(queue) > 0
        curr = queue[length(queue)]
        deleteat!(queue, length(queue))

        if curr.leaf push!(leaves, curr) end
        for c in getchildren(curr)
            push!(queue, c)
        end
    end

    return leaves
end

