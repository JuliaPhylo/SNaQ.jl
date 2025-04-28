# functions to undo incycle, containroot, gammaz
# originally in functions.jl
# Claudia March 2015


# ---------------------------------------- undo update of new hybridization --------------------------------


# function to undo updateInCycle which returns an array
# of edges/nodes changed
function undoInCycle!(edges::Array{Edge,1},nodes::Array{Node,1})
    for e in edges
        inCycle!(e, -1);
    end
    for n in nodes
        inCycle!(n, -1);
    end
end

# function to undo updateContainRoot (which returns an array
# of edges changed) and gives value bool, by default true
function undoContainRoot!(edges::Array{Edge,1}, bool::Bool)
    for e in edges
        e.containroot = bool
    end
end

undoContainRoot!(edges::Vector{Edge}) = undoContainRoot!(edges,true)

# function to undo updateGammaz which returns an array
# of edges changed
# it only changes the status of istIdentifiable to true
function undoistIdentifiable!(edges::Array{Edge,1})
    for e in edges
        !istIdentifiable(e) ? istIdentifiable!(e, true) : istIdentifiable!(e, false);
    end
end


"""
    undoGammaz!(node, network)

Undo `updateGammaz!` for the 2 cases: bad diamond I,II.
`node` should be a hybrid node.
Set length to edges that were not identifiable and
change edges' `gammaz` attribute to -1.0.
Recalculate branch lengths in terms of `gammaz`.  
*warning*: needs to know `incycle` attributes
"""
function undoGammaz!(node::Node, net::HybridNetwork)
    node.hybrid || error("cannot undo gammaz if starting node is not hybrid")
    if(isBadDiamondI(node))
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        other_maj = getOtherNode(edge_maj,node);
        other_min = getOtherNode(edge_min,node);
        edgebla,tree_edge_incycle1,tree_edge = hybridEdges(other_min);
        edgebla,tree_edge_incycle2,tree_edge = hybridEdges(other_maj);
        gammaz(other_min) != -1 || error("bad diamond I in node $(node.number) but no gammaz updated correctly")
        setLength!(tree_edge_incycle1,-log(1-gammaz(other_min)))
        gammaz(other_maj) != -1 || error("bad diamond I in node $(node.number) but no gammaz updated correctly")
        setLength!(tree_edge_incycle2,-log(1-gammaz(other_maj)))
        if approxEq(gammaz(other_maj),0.0) && approxEq(gammaz(other_min),0.0)
            setgamma!(edge_maj,0.0, true) # gamma could be anything if both gammaz are 0.0, but will set to 0.0
            setLength!(edge_maj,0.0)
            setLength!(edge_min,0.0)
        else
            setgamma!(edge_maj,gammaz(other_maj) / (gammaz(other_maj)+gammaz(other_min)), true)
        end
        gammaz!(other_min, -1.0)
        gammaz!(other_maj, -1.0)
        istIdentifiable!(tree_edge_incycle1, true);
        istIdentifiable!(tree_edge_incycle2, true);
        istIdentifiable!(edge_maj, true);
        istIdentifiable!(edge_min, true);
        isBadDiamondI!(node, false)
        numBad!(net, numBad(net) - 1)
    elseif(isBadDiamondII(node))
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        istIdentifiable!(tree_edge2, true)
        isBadDiamondII!(node, false)
    elseif(isBadTriangle(node))
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        istIdentifiable!(tree_edge2, true)
        isBadTriangle!(node, false)
    elseif(isVeryBadTriangle(node) || isExtBadTriangle(node))
        isVeryBadTriangle!(node, false)
        isExtBadTriangle!(node, false)
        hasVeryBadTriangle!(net, false)
    else
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        istIdentifiable!(edge_maj, isEdgeIdentifiable(edge_maj))
        istIdentifiable!(edge_min, isEdgeIdentifiable(edge_min))
    end
end
