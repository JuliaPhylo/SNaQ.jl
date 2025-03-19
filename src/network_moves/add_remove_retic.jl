using PhyloNetworks
import PhyloNetworks: addhybridedge!, directionalconflict, deletehybridedge!
import StatsBase: sample

const Edge = PhyloNetworks.Edge;
const Node = PhyloNetworks.Node;



"""
Splits `from_edge` and `to_edge` to create a new hybrid. Returns the new hybrid node.
"""
function add_hybrid!(from_edge::Edge, to_edge::Edge, N::HybridNetwork)
    is_valid_add_hybrid(from_edge, to_edge, N) || error("Invalid parameters given to `add_hybrid!`")
    newH, _ = addhybridedge!(N, from_edge, to_edge, true, 0.0, 0.25)
    getparent(from_edge).name = "i$(abs(rand(Int)) % 10000 + 1000)"
    newH.name = "H$(abs(rand(Int)) % 10000 + 1000)"
    return newH
end

function remove_hybrid!(hyb_node::Node, N::HybridNetwork, minor::Bool=true)
    rm_edge = minor ? getparentedgeminor(hyb_node) : getparentedge(hyb_node)
    deletehybridedge!(N, rm_edge)
end


"""
Randomly selects pairs of edges to use with `add_hybrid!` until a valid pair is found.
"""
function add_random_hybrid!(N::HybridNetwork, rng::TaskLocalRNG)

    niter::Int = 0
    e1::Edge = N.edge[1]    # placeholders
    e2::Edge = N.edge[1]

    while niter < 1000
        e1, e2 = sample(rng, N.edge, 2, replace=false)
        is_valid_add_hybrid(e1, e2, N) && return add_hybrid!(e1, e2, N)
        is_valid_add_hybrid(e2, e1, N) && return add_hybrid!(e2, e1, N)
    end

    error("Could not find any valid edge pairings to return.")

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







