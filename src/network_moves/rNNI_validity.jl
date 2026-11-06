
"""
Checks whether rNNI move of type `type` ([see Figure 4 here](https://doi.org/10.1371/journal.pcbi.1005611))
is valid with selected nodes `s,t,u,v`.
"""
is_valid_rNNI(s::Node, t::Node, u::Node, v::Node, type::Int) = type == 1 ? isvalidrNNI1(s, t, u, v) :
    type == 2 ? isvalidrNNI2(s, t, u, v) :
    type == 3 ? isvalidrNNI3(s, t, u, v) :
    type == 4 ? isvalidrNNI4(s, t, u, v) : error("`type` must be in [1, 2, 3, 4], received $(type).")

function areconnected(n1::Node, n2::Node)::Bool
    return any(e -> length(e.node) == 2 && (
        (e.node[1] == n1 && e.node[2] == n2) ||
        (e.node[1] == n2 && e.node[2] == n1)
    ), n1.edge)
end

function isvalidrNNI1(s::Node, t::Node, u::Node, v::Node)
    (s == t || s == u || s == v || t == u || t == v || u == v) && return false

    # If any of s, t, u, v are hybrids, then we need to enforce directionality
    if s.hybrid || t.hybrid || u.hybrid || v.hybrid
        cu = getchildren(u)
        (s in cu && v in cu) || return false

        !isdescendantof(s, v) || return false
        cv = getchildren(v)

        (!(s in cv) && !(v in getchildren(s))) || return false
    end

    # Directionality doesn't matter b/c we aren't actually rooted;
    # what matters is that the edges that need to be present are
    # indeed present (su, uv, tv) AND v is not a hybrid
    return areconnected(u, v) &&
        areconnected(s, u) &&
        areconnected(t, v)
end

function isvalidrNNI2(s::Node, t::Node, u::Node, v::Node)
    (s == t || s == u || s == v || t == u || t == v || u == v) && return false
    v in getchildren(u) || return false
    v in getchildren(t) || return false
    u in getchildren(s) || return false
    v.hybrid || return false
    !isdescendantof(t, u) || return false
    return true
end

function isvalidrNNI3(s::Node, t::Node, u::Node, v::Node)
    (s == t || s == u || s == v || t == u || t == v || u == v) && return false
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
function isvalidrNNI4(s::Node, t::Node, u::Node, v::Node)
    (s == t || s == u || s == v || t == u || t == v || u == v) && return false
    v in getchildren(u) || return false
    v in getchildren(t) || return false
    s in getchildren(u) || return false
    length(getchildren(u)) >= 2 || return false
    v.hybrid || return false
    !isdescendantof(u, t) || return false
    !isdescendantof(t, s) || return false
    return true
end


"""
Gathers all sets of nodes `s,t,u,v` (labels from [Fig 4 of this paper](https://doi.org/10.1371/journal.pcbi.1005611)) where `s,t,u,v` define a valid NNI move of type `type` (also defined in Fig 4).
"""
function allvalidrNNInodes(N::HybridNetwork, type::Int)
    if type == 1
        return allvalidrNNI1nodes(N)
    elseif type == 2
        return allvalidrNNI2nodes(N)
    elseif type == 3
        return allvalidrNNI3nodes(N)
    else
        return allvalidrNNI4nodes(N)
    end
end

function allvalidrNNI1nodes(N::HybridNetwork)::Vector{NTuple{4, Node}}
    nonleaves = [node for node in N.node if !node.leaf && !node.hybrid]
    stuv_combos = Vector{NTuple{4, Node}}();

    for i = 1:(length(nonleaves)-1)
        u = nonleaves[i]
        u.hybrid && continue
        for j = (i+1):length(nonleaves)
            v = nonleaves[j]
            v.hybrid && continue
            areconnected(u, v) || continue

            for se in u.edge
                s = getOtherNode(se, u)
                s == v && continue
                s.hybrid && continue
                for te in v.edge
                    t = getOtherNode(te, v)
                    t == u && continue
                    t.hybrid && continue
                    push!(stuv_combos, (s, t, u, v))
                end
            end
        end
    end
    return stuv_combos

    ####### Out-dated #########
    # valid_us = [node for node in N.node if !node.leaf && length(getchildren(node)) == 2];
    # length(valid_us) > 0 || return []

    # stuv_combos = [];
    # for u in valid_us
    #     u == getroot(N) && continue
    #     children = getchildren(u)
    #     (isdescendantof(children[1], children[2]) || isdescendantof(children[2], children[1])) && continue
    #     children[1] == children[2] && continue

    #     if !children[1].leaf && !(children[2] in getchildren(children[1]))
    #         for t in getchildren(children[1])
    #             children[2].hybrid && children[1] in getparents(children[2]) && continue
    #             push!(stuv_combos, (children[2], t, u, children[1]))
    #         end
    #     end

    #     if !children[2].leaf && !(children[1] in getchildren(children[2]))
    #         for t in getchildren(children[2])
    #             children[1].hybrid && children[2] in getparents(children[1]) && continue
    #             push!(stuv_combos, (children[1], t, u, children[2]))
    #         end
    #     end
    # end
    # all([isvalidrNNI1(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")

    # return stuv_combos
end

function allvalidrNNI2nodes(N::HybridNetwork)
    valid_vs = [node for node in N.node if node.hybrid];
    length(valid_vs) > 0 || return []

    stuv_combos = [];
    for v in valid_vs
        parents = getparents(v)
        parents[1] == parents[2] && continue
        p1_parents = getparents(parents[1])
        p2_parents = getparents(parents[2])
        
        if length(p1_parents) >= 1 && !isdescendantof(parents[2], parents[1])
            for gpa in p1_parents
                !isdescendantof(gpa, v) || continue
                (gpa == parents[2] || parents[1] == v) && continue
                push!(stuv_combos, (gpa, parents[2], parents[1], v))
            end
        end

        if length(p2_parents) >= 1 && !isdescendantof(parents[1], parents[2])
            for gpa in p2_parents
                !isdescendantof(gpa, v) || continue
                (gpa == parents[1] || parents[2] == v) && continue
                push!(stuv_combos, (gpa, parents[1], parents[2], v))
            end
        end
    end
    all([isvalidrNNI2(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")

    return stuv_combos
end

function allvalidrNNI3nodes(N::HybridNetwork)
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
    all([isvalidrNNI3(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")

    return stuv_combos
end

function allvalidrNNI4nodes(N::HybridNetwork)
    valid_vs = [node for node in N.node if node.hybrid];
    length(valid_vs) > 0 || return []

    stuv_combos = [];
    for v in valid_vs
        parents = getparents(v)
        p1_children = getchildren(parents[1])
        p2_children = getchildren(parents[2])

        if length(p1_children) == 2 && !isdescendantof(parents[1], parents[2])
            for s in p1_children
                s != parents[2] || continue
                s != v || continue
                isdescendantof(parents[2], s) && continue
                push!(stuv_combos, (s, parents[2], parents[1], v))
            end
        end

        if length(p2_children) == 2 && !isdescendantof(parents[2], parents[1])
            for s in p2_children
                s != parents[1] || continue
                s != v || continue
                isdescendantof(parents[1], s) && continue
                push!(stuv_combos, (s, parents[1], parents[2], v))
            end
        end
    end
    all([isvalidrNNI4(s, t, u, v) for (s, t, u, v) in stuv_combos]) || error("Some generated combos weren't valid...")
    return stuv_combos
end