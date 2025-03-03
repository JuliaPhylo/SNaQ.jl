using PhyloNetworks
import PhyloNetworks: getOtherNode


"""
Moves the origin of `hybrid` - equivalent to an rSPR move.
"""
function move_reticulate_origin!(hybrid::Node, new_origin::Edge, N::HybridNetwork)
    hybrid.hybrid || error("`hybrid` is not a hybrid node...")
    @debug "move_reticulate_origin!: ($(hybrid.name), $(getparent(new_origin).name) --> $(getchild(new_origin).name))"
    
    w = hybrid
    z = getparent(getparentedgeminor(hybrid))
    x, y = nothing, nothing

    if length(getparents(z)) > 0
        y = getchildren(z)[1] == hybrid ? getchildren(z)[2] : getchildren(z)[1]
        x = getparent(z)
    else
        # z is the "root" at the moment - but we ignore this notion
        # because the network is semi-directed
        x, y = [getOtherNode(edge, z) for edge in z.edge if getOtherNode(edge, z) != hybrid]
    end 

    xprime = getparent(new_origin)
    yprime = getchild(new_origin)

    is_valid_rSPR(w, x, y, z, xprime, yprime) || error("Invalid `move_reticulate_origin!` move.")
    return perform_rSPR!(N, w, x, y, z, xprime, yprime)
end

function move_random_reticulate_origin!(N::HybridNetwork)
    random_H_order = sample(1:N.numhybrids, N.numhybrids, replace=false)

    for H_idx in random_H_order
        H = N.hybrid[H_idx]
        random_edge_order = sample(1:N.numedges, N.numedges, replace=false)
        for edge_idx in random_edge_order
            is_valid_move_reticulate_origin(H, N.edge[edge_idx], N) || continue
            return move_reticulate_origin!(H, N.edge[edge_idx], N)
        end
    end

    error("No valid `move_random_reticulate_origin!` moves.")
end

function is_valid_move_reticulate_origin(hybrid::Node, new_origin::Edge, N::HybridNetwork)
    hybrid.hybrid || return false
    z = getparent(getparentedgeminor(hybrid))
    length(z.edge) == 3 || error("Node with only $(length(z.edge)) attached edges...")
    xprime = getparent(new_origin)
    yprime = getchild(new_origin)

    if length(getparents(z)) > 0
        y = getchildren(z)[1] == hybrid ? getchildren(z)[2] : getchildren(z)[1]
        x = getparent(z)
        return is_valid_rSPR(hybrid, x, y, z, xprime, yprime)
    else
        # z is the "root" at the moment - but we ignore this notion
        # because the network is semi-directed
        x, y = [getOtherNode(edge, z) for edge in z.edge if getOtherNode(edge, z) != hybrid]
        return is_valid_rSPR(hybrid, x, y, z, xprime, yprime)
    end 
end


function move_reticulate_target!(hybrid::Node, new_target::Edge, N::HybridNetwork)
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

function move_random_reticulate_target!(N::HybridNetwork)
    random_H_order = sample(1:N.numhybrids, N.numhybrids, replace=false)

    for H_idx in random_H_order
        H = N.hybrid[H_idx]
        random_edge_order = sample(1:N.numedges, N.numedges, replace=false)
        for edge_idx in random_edge_order
            is_valid_move_reticulate_target(H, N.edge[edge_idx], N) || continue
            return move_reticulate_target!(H, N.edge[edge_idx], N)
        end
    end

    error("No valid `move_random_reticulate_target!` moves.")
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