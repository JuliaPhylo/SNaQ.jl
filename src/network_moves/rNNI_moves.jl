using PhyloNetworks, StatsBase, Random
const Node = PhyloNetworks.Node;
const Edge = PhyloNetworks.Edge;
include("misc.jl")
include("rNNI_validity.jl")
include("../gradient_optimization/CF_recursive_blocks.jl")


"""
Performs an rNNI move on nodes `s,t,u,v` corresponding to those in Figure 4 of https://doi.org/10.1371/journal.pcbi.1005611.
Only implementing (1)-(4) WITHOUT (*) versions for now.
"""
function perform_rNNI!(N::HybridNetwork, s::Node, t::Node, u::Node, v::Node, type::Int)
    if type == 1
        return perform_rNNI1!(N, s, t, u, v)
    elseif type == 2
        return perform_rNNI2!(N, s, t, u, v)
    elseif type == 3
        return perform_rNNI3!(N, s, t, u, v)
    elseif type == 4
        return perform_rNNI4!(N, s, t, u, v)
    end
end


function perform_rNNI1!(N::HybridNetwork, s::Node, t::Node, u::Node, v::Node)
    is_valid_rNNI1(s, t, u, v) || error("Topological conditions for rNNI(1) not met.")
    @debug "MOVE: rNNI(1) - $((s.name, t.name, u.name, v.name))"

    # u: loses s as a child and gains t as a child
    # v: loses t as a child and gains s as a child
    su_idx::Int = findfirst(e -> u in e.node, s.edge)
    tv_idx::Int = findfirst(e -> v in e.node, t.edge)
    edge_su::Edge = s.edge[su_idx]
    edge_tv::Edge = t.edge[tv_idx]
    replace!(edge_su.node, s => t)
    replace!(edge_tv.node, t => s)
    replace!(s.edge, edge_su => edge_tv)
    replace!(t.edge, edge_tv => edge_su)
    #edge_sv::Edge = s.edge[findfirst(e -> v in e.node, s.edge)]

    # temp_number = edge_su.number
    # edge_su.number = edge_tv.number
    # edge_tv.number = temp_number

    # Swap the hybridization-related info of these edges. That way, if either
    # `s` or `t` are hybrids, they maintain their status
    swap_edge_hybrid_info!(edge_su, edge_tv)
    return 1
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
    return 2
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
    v.name = v.name == "" ? "H$(abs(rand(Int) % 10000) + 100)" : v.name[1] == "H" ? v.name : "H$(v.name)"
    replace!(N.hybrid, u => v)
    return 3
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
    u.name = u.name == "" ? "H$(abs(rand(Int) % 10000) + 100)" : u.name[1] == "H" ? u.name : "H$(u.name)"
    replace!(N.hybrid, v => u)
    return 4
end


"""
Samples nodes `s,t,u,v` from `N` to perform the rNNI(X) move where X is given
by the variable `type`.
"""
function sample_rNNI_parameters(N::HybridNetwork, type::Int, rng::TaskLocalRNG)
    1 ≤ type ≤ 4 || error("`type` must be 1, 2, 3, or 4.")

    valid_stuvs = type == 1 ? all_valid_rNNI1_nodes(N) :
        type == 2 ? all_valid_rNNI2_nodes(N) :
        type == 3 ? all_valid_rNNI3_nodes(N) :
        all_valid_rNNI4_nodes(N)
    
    length(valid_stuvs) == 0 && return nothing
    return sample(rng, valid_stuvs)
end


"""
Add a docstring later - I'm going on a walk 
"""
function sample_random_rNNI(N::HybridNetwork, rng::TaskLocalRNG, probs::Vector{<:Real}=[0.25, 0.25, 0.25, 0.25])
    (length(probs) == 4 && sum(probs) ≈ 1) || error("`probs` must have length 4 and sum to 1.")
    
    r = rand()
    if r < probs[1]
        valid_stuvs = all_valid_rNNI1_nodes(N)
        if length(valid_stuvs) == 0
            p234 = probs[2] + probs[3] + probs[4]
            return perform_random_rNNI!(N, rng, [0, probs[2], probs[3], probs[4]] ./ p234)
        end
        s, t, u, v = sample(rng, valid_stuvs)
        return perform_rNNI1!(N, s, t, u, v)
    elseif r < probs[1] + probs[2]
        valid_stuvs = all_valid_rNNI2_nodes(N)
        if length(valid_stuvs) == 0
            p134 = probs[1] + probs[3] + probs[4]
            return perform_random_rNNI!(N, rng, [probs[1], 0, probs[3], probs[4]] ./ p134)
        end
        s, t, u, v = sample(rng, valid_stuvs)
        return perform_rNNI2!(N, s, t, u, v)
    elseif r < probs[1] + probs[2] + probs[3]
        valid_stuvs = all_valid_rNNI3_nodes(N)
        if length(valid_stuvs) == 0
            p124 = probs[1] + probs[2] + probs[4]
            return perform_random_rNNI!(N, rng, [probs[1], probs[2], 0, probs[4]] ./ p124)
        end
        s, t, u, v = sample(rng, valid_stuvs)
        return perform_rNNI3!(N, s, t, u, v)
    else
        valid_stuvs = all_valid_rNNI4_nodes(N)
        if length(valid_stuvs) == 0
            p123 = probs[1] + probs[2] + probs[3]
            return perform_random_rNNI!(N, rng, [probs[1], probs[2], probs[3], 0] ./ p123)
        end
        s, t, u, v = sample(rng, valid_stuvs)
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