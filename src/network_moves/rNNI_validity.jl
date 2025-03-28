using PhyloNetworks, StatsBase
const Node = PhyloNetworks.Node;
const Edge = PhyloNetworks.Edge;
include("misc.jl")

"""
Checks whether rNNI move of type `type` ([see Figure 4 here](https://doi.org/10.1371/journal.pcbi.1005611))
is valid with selected nodes `s,t,u,v`.
"""
is_valid_rNNI(s::Node, t::Node, u::Node, v::Node, type::Int) = type == 1 ? is_valid_rNNI1(s, t, u, v) :
    type == 2 ? is_valid_rNNI2(s, t, u, v) :
    type == 3 ? is_valid_rNNI3(s, t, u, v) :
    type == 4 ? is_valid_rNNI4(s, t, u, v) : error("`type` must be in [1, 2, 3, 4], received $(type).")


function is_valid_rNNI1(s::Node, t::Node, u::Node, v::Node)
    !v.hybrid || return false

    cu = getchildren(u)
    (s in cu && v in cu) || return false

    cv = getchildren(v)
    t in cv || return false

    (!(s in cv) && !(v in getchildren(s))) || return false
    return true
end

function is_valid_rNNI2(s::Node, t::Node, u::Node, v::Node)
    v in getchildren(u) || return false
    v in getchildren(t) || return false
    u in getchildren(s) || return false
    v.hybrid || return false
    !is_descendant_of(t, u) || return false
    return true
end

function is_valid_rNNI3(s::Node, t::Node, u::Node, v::Node)
    !v.hybrid || return false
    t in getchildren(v) || return false
    u in getchildren(s) || return false
    v in getchildren(u) || return false
    u.hybrid || return false
    getparents(u)[1] != getparents(u)[2] || return false
    return true
end

# Basic visual relationships in Figure 4, but ALSO
# `t` must NOT be a descendant of `u`
function is_valid_rNNI4(s::Node, t::Node, u::Node, v::Node)
    v in getchildren(u) || return false
    v in getchildren(t) || return false
    s in getchildren(u) || return false
    length(getchildren(u)) == 2 || return false
    v.hybrid || return false
    !is_descendant_of(u, t) || return false
    !is_descendant_of(t, s) || return false
    return true
end


"""
Gathers all sets of nodes `s,t,u,v` (labels from [Fig 4 of this paper](https://doi.org/10.1371/journal.pcbi.1005611)) where `s,t,u,v` define a valid NNI move of type `type` (also defined in Fig 4).
"""
function all_valid_rNNI_nodes(N::HybridNetwork, type::Int)
    if type == 1
        return all_valid_rNNI1_nodes(N)
    elseif type == 2
        return all_valid_rNNI2_nodes(N)
    elseif type == 3
        return all_valid_rNNI3_nodes(N)
    else
        return all_valid_rNNI4_nodes(N)
    end
end

function sample_valid_rNNI_nodes(N::HybridNetwork; probs::Vector{<:Real}=[0.25, 0.25, 0.25, 0.25])
    (length(probs) == 4 && sum(probs) ≈ 1) || error("`probs` must have length 4 and sum to 1.")
    
    r = rand()
    if r < probs[1]
        return sample_valid_rNNI1_nodes(N)
    elseif r < probs[1] + probs[2]
        return sample_valid_rNNI2_nodes(N)
    elseif r < probs[1] + probs[2] + probs[3]
        return sample_valid_rNNI3_nodes(N)
    else
        return sample_valid_rNNI4_nodes(N)
    end
end

function all_valid_rNNI1_nodes(N::HybridNetwork)
    valid_us = [node for node in N.node if !node.leaf && length(getchildren(node)) == 2];
    length(valid_us) > 0 || return []

    stuv_combos = [];
    for u in valid_us
        u == getroot(N) && continue
        children = getchildren(u)

        if !children[1].leaf && !(children[2] in getchildren(children[1])) && !children[1].hybrid
            for t in getchildren(children[1])
                children[2].hybrid && children[1] in getparents(children[2]) && continue
                push!(stuv_combos, (children[2], t, u, children[1]))
            end
        end

        if !children[2].leaf && !(children[1] in getchildren(children[2])) && !children[2].hybrid
            for t in getchildren(children[2])
                children[1].hybrid && children[2] in getparents(children[1]) && continue
                push!(stuv_combos, (children[1], t, u, children[2]))
            end
        end
    end
    all([is_valid_rNNI1(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")

    return stuv_combos
end

function all_valid_rNNI2_nodes(N::HybridNetwork)
    valid_vs = [node for node in N.node if node.hybrid];
    length(valid_vs) > 0 || return []

    stuv_combos = [];
    for v in valid_vs
        parents = getparents(v)
        p1_parents = getparents(parents[1])
        p2_parents = getparents(parents[2])
        
        if length(p1_parents) >= 1 && !is_descendant_of(parents[2], parents[1])
            for gpa in p1_parents
                !is_descendant_of(gpa, v) || continue
                push!(stuv_combos, (gpa, parents[2], parents[1], v))
            end
        end

        if length(p2_parents) >= 1 && !is_descendant_of(parents[1], parents[2])
            for gpa in p2_parents
                !is_descendant_of(gpa, v) || continue
                push!(stuv_combos, (gpa, parents[1], parents[2], v))
            end
        end
    end
    all([is_valid_rNNI2(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")

    return stuv_combos
end

function all_valid_rNNI3_nodes(N::HybridNetwork)
    valid_us = [node for node in N.node if node.hybrid && !getchild(node).leaf];
    length(valid_us) > 0 || return []

    stuv_combos = [];
    for u in valid_us
        v = getchild(u)
        v.hybrid && continue
        ss = getparents(u)
        ss[1] != ss[2] || continue
        
        for s in ss
            for t in getchildren(v)
                push!(stuv_combos, (s, t, u, v))
            end
        end
    end
    all([is_valid_rNNI3(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")

    return stuv_combos
end

function all_valid_rNNI4_nodes(N::HybridNetwork)
    valid_vs = [node for node in N.node if node.hybrid];
    length(valid_vs) > 0 || return []

    stuv_combos = [];
    for v in valid_vs
        parents = getparents(v)
        p1_children = getchildren(parents[1])
        p2_children = getchildren(parents[2])

        if length(p1_children) == 2 && !is_descendant_of(parents[1], parents[2])
            for s in p1_children
                s != parents[2] || continue
                s != v || continue
                is_descendant_of(parents[2], s) && continue
                push!(stuv_combos, (s, parents[2], parents[1], v))
            end
        end

        if length(p2_children) == 2 && !is_descendant_of(parents[2], parents[1])
            for s in p2_children
                s != parents[1] || continue
                s != v || continue
                is_descendant_of(parents[1], s) && continue
                push!(stuv_combos, (s, parents[1], parents[2], v))
            end
        end
    end
    all([is_valid_rNNI4(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")
    return stuv_combos
end