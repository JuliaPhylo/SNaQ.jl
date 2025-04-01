using PhyloNetworks
import PhyloNetworks: addhybridedge!, directionalconflict, deletehybridedge!
import StatsBase: sample

const Edge = PhyloNetworks.Edge;
const Node = PhyloNetworks.Node;



"""
Splits `from_edge` and `to_edge` to create a new hybrid. Returns the new hybrid node.
"""
function add_hybrid!(N::HybridNetwork, from_edge::Edge, to_edge::Edge)
    is_valid_add_hybrid(from_edge, to_edge, N) || error("Invalid parameters given to `add_hybrid!`")
    newH, _ = addhybridedge!(N, from_edge, to_edge, true, 0.0, 0.25)
    getparent(from_edge).name = "i$(abs(rand(Int)) % 100000 + 1000)"
    newH.name = "H$(abs(rand(Int)) % 100000 + 1000)"
    return newH
end


"""
Removes hybrid node `hyb_node` from network `N`. Removes the minor edge if
argument `minor` is `true`, else removes the `major` edge.
"""
function remove_hybrid!(N::HybridNetwork, hyb_node::Node, minor::Bool=true)
    rm_edge = minor ? getparentedgeminor(hyb_node) : getparentedge(hyb_node)
    deletehybridedge!(N, rm_edge)
end


"""
Randomly samples a hybrid node from network `N`.
"""
function sample_remove_hybrid_parameters(N::HybridNetwork, rng::TaskLocalRNG)
    return sample(rng, N.hybrid)
end


"""
Randomly selects pairs of edges to use with `add_hybrid!` until a valid pair is found.
Returns that pair of nodes, or `nothing` if no valid pair is found (this should
never happen).
"""
function sample_add_hybrid_parameters(N::HybridNetwork, rng::TaskLocalRNG)

    niter::Int = 0
    e1::Edge = N.edge[1]    # placeholders
    e2::Edge = N.edge[1]

    while niter < 1000
        e1, e2 = sample(rng, N.edge, 2, replace=false)
        is_valid_add_hybrid(e1, e2, N) && return (e1, e2)
        is_valid_add_hybrid(e2, e1, N) && return (e2, e1)
    end

    return nothing

end


"""
Checks whether the `add_hybrid!` move defined by this `from_edge` and
`to_edge` choice is valid. A move must satisfy the following:

1. must not create a direction conflict
"""
function is_valid_add_hybrid(from_edge::Edge, to_edge::Edge, N::HybridNetwork)::Bool
    # Condition (1)
    !directionalconflict(getparent(from_edge), to_edge, true) || return false

    return true
end


# INVALID example







