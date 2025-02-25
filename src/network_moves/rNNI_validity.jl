using PhyloNetworks
const Node = PhyloNetworks.Node;
const Edge = PhyloNetworks.Edge;

"""
Checks whether rNNI move of type `type` ([see Figure 4 here](https://doi.org/10.1371/journal.pcbi.1005611))
is valid with selected nodes `s,t,u,v`.
"""
is_valid_rNNI(s::Node, t::Node, u::Node, v::Node, type::Int) = type == 1 ? is_valid_rNNI1(s, t, u, v) :
    type == 2 ? is_valid_rNNI2(s, t, u, v) :
    type == 3 ? is_valid_rNNI3(s, t, u, v) :
    type == 4 ? is_valid_rNNI4(s, t, u, v) : error("`type` must be in [1, 2, 3, 4], received $(type).")



function is_valid_rNNI1(s::Node, t::Node, u::Node, v::Node)
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
    return true
end

function is_valid_rNNI3(s::Node, t::Node, u::Node, v::Node)
    t in getchildren(v) || return false
    u in getchildren(s) || return false
    v in getchildren(u) || return false
    u.hybrid || return false
    return true
end

function is_valid_rNNI4(s::Node, t::Node, u::Node, v::Node)
    v in getchildren(u) || return false
    v in getchildren(t) || return false
    s in getchildren(u) || return false
    v.hybrid || return false
    return true
end