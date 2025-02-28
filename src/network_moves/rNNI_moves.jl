using PhyloNetworks, StatsBase
const Node = PhyloNetworks.Node;
const Edge = PhyloNetworks.Edge;
include("misc.jl")
include("rNNI_validity.jl")


"""
Performs an rNNI move on nodes `s,t,u,v` corresponding to those in Figure 4 of https://doi.org/10.1371/journal.pcbi.1005611.
Only implementing (1)-(4) WITHOUT (*) versions for now.
"""
function perform_rNNI!(N::HybridNetwork, s::Node, t::Node, u::Node, v::Node, type::Int)
    if type == 1
        perform_rNNI1!(N, s, t, u, v)
    elseif type == 2
        perform_rNNI2!(N, s, t, u, v)
    elseif type == 3
        perform_rNNI3!(N, s, t, u, v)
    elseif type == 4
        perform_rNNI4!(N, s, t, u, v)
    end
end


function perform_rNNI1!(N::HybridNetwork, s::Node, t::Node, u::Node, v::Node)
    is_valid_rNNI1(s, t, u, v) || error("Topological conditions for rNNI(1) not met.")
    @debug "MOVE: rNNI(1) - $((s.name, t.name, u.name, v.name))"

    # u: loses s as a child and gains t as a child
    # v: loses t as a child and gains s as a child
    edge_su::Edge = s.edge[findfirst(e -> u in e.node, s.edge)]
    edge_tv::Edge = t.edge[findfirst(e -> v in e.node, t.edge)]
    replace!(edge_su.node, s => t)
    replace!(edge_tv.node, t => s)
    replace!(s.edge, edge_su => edge_tv)
    replace!(t.edge, edge_tv => edge_su)

    # Swap the hybridization-related info of these edges. That way, if either
    # `s` or `t` are hybrids, they maintain their status
    swap_edge_hybrid_info!(edge_su, edge_tv)
end


function perform_rNNI2!(N::HybridNetwork, s::Node, t::Node, u::Node, v::Node)
    is_valid_rNNI2(s, t, u, v) || error("Topological conditions for rNNI(2) not met.")
    @debug "MOVE: rNNI(2) - $((s.name, t.name, u.name, v.name))"

    # s: loses u as a child and gains v as a child
    # t: loses v as a child and gains u as a child
    edge_su::Edge = s.edge[findfirst(e -> u in e.node, s.edge)]
    edge_tv::Edge = t.edge[findfirst(e -> v in e.node, t.edge)]
    edge_uv::Edge = u.edge[findfirst(e -> v in e.node, u.edge)]
    (edge_uv.hybrid && edge_tv.hybrid) || error("Topological conditions for move not met")

    replace!(edge_su.node, u => v)
    replace!(edge_tv.node, v => u)
    replace!(u.edge, edge_su => edge_tv)
    replace!(v.edge, edge_tv => edge_su)

    # Swap the hybridization-related info of these edges. That way, if
    # `u` is a hybrid, it maintains its status, and `v` maintains its status as a hybrid
    swap_edge_hybrid_info!(edge_su, edge_tv)
end


function perform_rNNI3!(N::HybridNetwork, s::Node, t::Node, u::Node, v::Node)
    is_valid_rNNI3(s, t, u, v) || error("Topological conditions for rNNI(3) not met.")
    @debug "MOVE: rNNI(3) - $((s.name, t.name, u.name, v.name))"

    # s: loses u as a child and gains v as a child
    # t: loses v as a parent and gains u as a parent
    edge_su::Edge = s.edge[findfirst(e -> u in e.node, s.edge)]
    edge_tv::Edge = t.edge[findfirst(e -> v in e.node, t.edge)]
    edge_uv::Edge = u.edge[findfirst(e -> v in e.node, u.edge)]
    edge_u::Edge = u.edge[findfirst(e -> !(s in e.node) && !(v in e.node), u.edge)]
    edge_su.hybrid || error("Topological conditions for move not met")

    replace!(edge_su.node, u => v)
    replace!(edge_tv.node, v => u)
    replace!(u.edge, edge_su => edge_tv)
    replace!(v.edge, edge_tv => edge_su)

    # New reticulate relationships:
    # - `u` no longer a hybrid
    # - `v` now a hybrid
    # - `s,t` should be unchanged
    #   to any of `s,t,v` is now a hybrid edge

    # This is the edge in Figure 4 (3) that is parent to `u` and black (not gray)
    swap_edge_hybrid_info!(edge_u, edge_uv)
    u.hybrid = false
    v.hybrid = true
    replace!(N.hybrid, u => v)
end


function perform_rNNI4!(N::HybridNetwork, s::Node, t::Node, u::Node, v::Node)
    is_valid_rNNI4(s, t, u, v) || error("Topological conditions for rNNI(4) not met.")
    @debug "MOVE: rNNI(4) - $((s.name, t.name, u.name, v.name))"

    # t: loses v as a child and gains u as a child
    # s: loses u as a parent and gains v as a parent
    edge_su::Edge = s.edge[findfirst(e -> u in e.node, s.edge)]
    edge_tv::Edge = t.edge[findfirst(e -> v in e.node, t.edge)]
    edge_uv::Edge = u.edge[findfirst(e -> v in e.node, u.edge)]
    edge_u::Edge = u.edge[findfirst(e -> !(s in e.node) && !(v in e.node), u.edge)]

    replace!(edge_su.node, u => v)
    replace!(edge_tv.node, v => u)
    replace!(u.edge, edge_su => edge_tv)
    replace!(v.edge, edge_tv => edge_su)

    # New reticulate relationships:
    # - `v` no longer a hybrid
    # - `u` now a hybrid
    # - `s,t` should be unchanged
    # - the edge parent to `u` that is not connected
    #   to any of `s,t,v` is now a hybrid edge
    swap_edge_hybrid_info!(edge_u, edge_uv)
    u.hybrid = true
    v.hybrid = false
    replace!(N.hybrid, v => u)
end


"""
Add a docstring later - I'm going on a walk 
"""
function perform_random_rNNI!(N::HybridNetwork, probs::Vector{<:Real}=[0.9, 0.1/3, 0.1/3, 0.1/3])
    (length(probs) == 4 && sum(probs) ≈ 1) || error("`probs` must have length 4 and sum to 1.")
    
    r = rand()
    if r < probs[1]
        valid_stuvs = all_valid_rNNI1_nodes(N)
        if length(valid_stuvs) == 0
            p234 = probs[2] + probs[3] + probs[4]
            return perform_random_rNNI!(N, [0, probs[2], probs[3], probs[4]] ./ p234)
        end
        s, t, u, v = sample(valid_stuvs)
        return perform_rNNI1!(N, s, t, u, v)
    elseif r < probs[1] + probs[2]
        valid_stuvs = all_valid_rNNI2_nodes(N)
        if length(valid_stuvs) == 0
            p134 = probs[1] + probs[3] + probs[4]
            return perform_random_rNNI!(N, [probs[1], 0, probs[3], probs[4]] ./ p134)
        end
        s, t, u, v = sample(valid_stuvs)
        return perform_rNNI2!(N, s, t, u, v)
    elseif r < probs[1] + probs[2] + probs[3]
        valid_stuvs = all_valid_rNNI3_nodes(N)
        if length(valid_stuvs) == 0
            p124 = probs[1] + probs[2] + probs[4]
            return perform_random_rNNI!(N, [probs[1], probs[2], 0, probs[4]] ./ p124)
        end
        s, t, u, v = sample(valid_stuvs)
        return perform_rNNI3!(N, s, t, u, v)
    else
        valid_stuvs = all_valid_rNNI4_nodes(N)
        if length(valid_stuvs) == 0
            p123 = probs[1] + probs[2] + probs[3]
            return perform_random_rNNI!(N, [probs[1], probs[2], probs[3], 0] ./ p123)
        end
        s, t, u, v = sample(valid_stuvs)
        return perform_rNNI4!(N, s, t, u, v)
    end
end


# EDGE attributes:
# - number:         unique ID
# - node:           attached nodes, always 2 in our case
# - ischild1:       TRUE if node[1] is child of the edge, FALSE if node[1] is parent of the edge
# - length:         branch length
# - hybrid:         whether the edge is a tree edge or hybrid edge
# - gamma:          γ
# - ismajor:        whether the edge is the major path to the child node
# - containroot:    is the interior of this edge a valid rooting position?





# NODE attributes:
# - edge:   attached edges (order does not matter)
# - hybrid: is the node a hybrid
# - leaf:   is the node a leaf
# - name:   node's name
# - number: unique ID
# 
# these ones can be ignored:
# - booln1-6
# - int8n3
# - intn1-2
# - fvalue
# - prev