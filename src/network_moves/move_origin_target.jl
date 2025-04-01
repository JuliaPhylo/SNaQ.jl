using PhyloNetworks
using StatsBase
import PhyloNetworks: getOtherNode


"""
Selects a random hybrid node and a random edge that is very near the origin of
that hybrid to perform the `move_reticulate_origin` move on. Returns
the selected parameters.
"""
function sample_move_reticulate_origin_local_parameters(N::HybridNetwork, rng::TaskLocalRNG)
    N.numhybrids > 0 || error("`N` must have at least 1 hybrid node.")
    
    # Cycle through each hybrid in a random order, that way if 1 hybrid
    # does not have any valid moves, another might.
    hyb_perm = randperm(rng, N.numhybrids)
    for hybrid in N.hybrid[hyb_perm]
        minor_edge = getparentedgeminor(hybrid)
        origin = getparent(minor_edge)

        # Edges near `hybrid` that might be valid for the move
        candidate_edges = reduce(vcat, [[next_edge for next_edge in getOtherNode(e, origin).edge if next_edge != e] for e in origin.edge if e != minor_edge])

        for j in sample(rng, 1:length(candidate_edges), length(candidate_edges), replace=false)
            is_valid_move_reticulate_origin(hybrid, candidate_edges[j], N) && return (hybrid, candidate_edges[j])
        end
    end

    return nothing
end


"""
Selects a random hybrid node and a random edge that is very near the origin of
that hybrid to perform the `move_reticulate_target` move on. Returns
the selected parameters.
"""
function sample_move_reticulate_target_local_parameters(N::HybridNetwork, rng::TaskLocalRNG)
    N.numhybrids > 0 || error("`N` must have at least 1 hybrid node.")

    # Cycle through each hybrid in a random order, that way if 1 hybrid
    # does not have any valid moves, another might.
    hyb_perm = randperm(rng, N.numhybrids)
    for hybrid in N.hybrid[hyb_perm]
        major_parent = getparent(getparentedge(hybrid))
        candidate_edges = [e for e in major_parent.edge if e != getparentedge(hybrid)]
        
        hyb_children = getchildren(hybrid)
        hyb_children = (length(hyb_children) == 1) ? getchildren(hyb_children[1]) : hyb_children
        for child in hyb_children
            push!(candidate_edges, getparentedge(child))
        end

        for j in sample(rng, 1:length(candidate_edges), length(candidate_edges), replace=false)
            is_valid_move_reticulate_target(hybrid, candidate_edges[j], N) && return (hybrid, candidate_edges[j])
        end
    end
    
    return nothing
end


"""
Moves the origin of `hybrid` - equivalent to an rSPR move.
"""
function move_reticulate_origin!(N::HybridNetwork, hybrid::Node, new_origin::Edge)
    hybrid.hybrid || error("`hybrid` is not a hybrid node...")
    @debug "move_reticulate_origin!: ($(hybrid.name), $(getparent(new_origin).name) --> $(getchild(new_origin).name))"
    
    w = hybrid
    z = getparent(getparentedgeminor(hybrid))
    x, y = nothing, nothing

    if length(getparents(z)) > 0
        y = length(getchildren(z)) == 1 || getchildren(z)[1] != hybrid ? getchildren(z)[1] : getchildren(z)[2]
        x = getparent(z)
    else
        # z is the "root" at the moment - but we ignore this notion
        # because the network is semi-directed
        other_nodes = [getOtherNode(edge, z) for edge in z.edge if getOtherNode(edge, z) != hybrid]
        length(other_nodes) == 2 || return false
        x, y = other_nodes
    end 

    xprime = getparent(new_origin)
    yprime = getchild(new_origin)

    is_valid_rSPR(w, x, y, z, xprime, yprime) || error("Invalid `move_reticulate_origin!` move.")
    return perform_rSPR!(N, w, x, y, z, xprime, yprime)
end


"""
Samples a random hybrid from `N` and a random edge in `N` from which to move
the origin of the selected hybrid to.
"""
function sample_move_reticulate_origin_parameters(N::HybridNetwork, rng::TaskLocalRNG)
    random_H_order = sample(rng, 1:N.numhybrids, N.numhybrids, replace=false)

    for H_idx in random_H_order
        H = N.hybrid[H_idx]
        random_edge_order = sample(rng, 1:N.numedges, N.numedges, replace=false)
        for edge_idx in random_edge_order
            is_valid_move_reticulate_origin(H, N.edge[edge_idx], N) || continue
            return (H, N.edge[edge_idx])
        end
    end

    return nothing
end

function is_valid_move_reticulate_origin(hybrid::Node, new_origin::Edge, N::HybridNetwork)
    hybrid.hybrid || return false
    z = getparent(getparentedgeminor(hybrid))
    length(z.edge) == 3 || error("Node with only $(length(z.edge)) attached edges...")
    xprime = getparent(new_origin)
    yprime = getchild(new_origin)

    if length(getparents(z)) > 0
        y = length(getchildren(z)) == 1 || getchildren(z)[1] != hybrid ? getchildren(z)[1] : getchildren(z)[2]
        x = getparent(z)
        return is_valid_rSPR(hybrid, x, y, z, xprime, yprime)
    else
        # z is the "root" at the moment - but we ignore this notion
        # because the network is semi-directed
        other_nodes = [getOtherNode(edge, z) for edge in z.edge if getOtherNode(edge, z) != hybrid]
        length(other_nodes) == 2 || return false
        x, y = other_nodes
        return is_valid_rSPR(hybrid, x, y, z, xprime, yprime)
    end 
end


function move_reticulate_target!(N::HybridNetwork, hybrid::Node, new_target::Edge)
    hybrid.hybrid || error("`hybrid` is not a hybrid node...")
    @debug "move_reticulate_target!: ($(hybrid.name), $(getparent(new_target).name) --> $(getchild(new_target).name))"

    x = getparent(getparentedge(hybrid))
    w = getparent(getparentedgeminor(hybrid))
    y = getchild(hybrid)

    xprime = getparent(new_target)
    yprime = getchild(new_target)

    is_valid_rSPR(w, x, y, hybrid, xprime, yprime) || error("Invalid `move_reticulate_target!` move.")
    return perform_rSPR!(N, w, x, y, hybrid, xprime, yprime)
end


"""
Samples a random hybrid from `N` and a random edge in `N` from which to move
the target of the selected hybrid to.
"""
function sample_move_reticulate_target_parameters(N::HybridNetwork, rng::TaskLocalRNG)
    random_H_order = sample(rng, 1:N.numhybrids, N.numhybrids, replace=false)

    for H_idx in random_H_order
        H = N.hybrid[H_idx]
        random_edge_order = sample(rng, 1:N.numedges, N.numedges, replace=false)
        for edge_idx in random_edge_order
            is_valid_move_reticulate_target(H, N.edge[edge_idx], N) || continue
            return (H, N.edge[edge_idx])
        end
    end

    return nothing
end

function is_valid_move_reticulate_target(hybrid::Node, new_target::Edge, N::HybridNetwork)
    hybrid.hybrid || error("`hybrid` is not a hybrid node...")
    x = getparent(getparentedge(hybrid))
    w = getparent(getparentedgeminor(hybrid))
    y = getchild(hybrid)

    xprime = getparent(new_target)
    yprime = getchild(new_target)

    return is_valid_rSPR(w, x, y, hybrid, xprime, yprime)
end