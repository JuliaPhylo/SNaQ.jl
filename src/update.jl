# functions to update incycle, containroot, gammaz
# Claudia March 2015
#####################

# --------------------------------------- update incycle, root, gammaz -------------------------------------------

# function to update inCycle (with priority queue) after becoming part of a network
# based on program 3 CS367 with priority queue
# expected to be much faster than the other two udpateInCycle (queue and recursive)
# input: hybrid node around which we want to update inCycle
# returns tuple: flag, nocycle, array of edges changed, array of nodes changed
#   flag: false if cycle intersects existing cycle or number of nodes in cycle < 3
#         (there is the possibility of returning edges in intersection: path)
#         true if cycle does not intersect existing cycle
#   nocycle: true if there is no cycle (instead of error). it is used in addHybridization
# calculates also the number of nodes in the cycle and put as hybrid node attribute "k"
# warning: it is not checking if hybrid node or minor hybrid edge
#          were already part of a cycle (inCycle = -1)
#          But it is checking so for the other edges in cycle
# warning: it needs extra things: visited attribute, prev attribute
#          unlike updateInCycle recursive, but it is expected
#          to be much faster
function updateInCycle!(net::HybridNetwork,node::Node)
    node.hybrid || error("node is not hybrid")
    start = node
    inCycle!(node, node.number)
    k!(node, 1)
    hybedge = getparentedgeminor(node)
    inCycle!(hybedge, node.number)
    lastnode = getOtherNode(hybedge,node)
    dist = 0
    queue = PriorityQueue()
    path = Node[]
    edges_changed!(net, Edge[])
    nodes_changed!(net, Node[])
    push!(edges_changed(net),hybedge)
    push!(nodes_changed(net),node)
    found = false
    visited!(net, falses(length(net.node)))
    enqueue!(queue,node,dist)
    while !found
        if isempty(queue)
            return false, true, edges_changed(net), nodes_changed(net)
        end
        curr = dequeue!(queue)
        if isEqual(curr,lastnode)
            found = true
            push!(path,curr)
        elseif !visited(net)[getIndex(curr,net)]
            visited(net)[getIndex(curr,net)] = true
            atstart = isEqual(curr,start)
            for e in curr.edge
                e.ismajor || continue
                other = getOtherNode(e,curr)
                if atstart || (!other.leaf && !visited(net)[getIndex(other,net)])
                    other.prev = curr
                    dist = dist+1
                    enqueue!(queue,other,dist)
                end
            end
        end
    end # end while
    curr = pop!(path)
    while !isEqual(curr, start)
        if inCycle(curr)!= -1
            push!(path,curr)
            curr = curr.prev
        else
            inCycle!(curr, start.number)
            push!(nodes_changed(net), curr)
            k!(node, k(node) + 1)
            edge = getconnectingedge(curr, curr.prev)
            inCycle!(edge, start.number)
            push!(edges_changed(net), edge)
            curr = curr.prev
        end
    end
    flag = isempty(path) # || k(node)<3
    !flag && @debug "warning: new cycle intersects existing cycle"
    return flag, false, edges_changed(net), nodes_changed(net)
end

"""
    updateContainRoot!(HybridNetwork, Node)
    traverseContainRoot!(Node, Edge, edges_changed::Array{Edge,1}, rightDir::Vector{Bool})

The input `node` to `updateContainRoot!` must be a hybrid node
(can come from PhyloNetworks.searchHybridNode).
`updateContainRoot!` starts at the input node and calls `traverseContainRoot!`,
which traverses the network recursively.
By default, containroot attributes of edges are true.
Changes `containroot` to false for all the visited edges: those
below the input node, but not beyond any other hybrid node.

`updateContainRoot!` Returns a `flag` and an array of edges whose
containroot has been changed from true to false.
`flag` is false if the set of edges to place the root is empty

In `traverseContainRoot!`, `rightDir` turns false if hybridizations
have incompatible directions (vector of length 1, to be modified).

Warning:

- does *not* update `containroot` of minor hybrid edges.
- assumes correct `ismajor` attributes: to stop the recursion at minor hybrid edges.
- assumes correct hybrid attributes of both nodes & edges: to check if various
  hybridizations have compatible directions.
  For each hybrid node that is encountered, checks if it was reached
  via a hybrid edge (ok) or tree edge (not ok).

`rightDir`: vector of length 1 boolean, to be mutable and modified by the function
"""
function traverseContainRoot!(node::Node, edge::Edge, edges_changed::Array{Edge,1}, rightDir::Vector{Bool})
    if node.hybrid
        if edge.hybrid
            edge.ismajor || error("hybrid edge $(edge.number) is minor and we should not traverse the graph through minor edges")
            DEBUGC && @debug "traverseContainRoot reaches hybrid node $(node.number) through major hybrid edge $(edge.number)"
            rightDir[1] &= true  # This line has no effect: x && true = x
        else #approach hybrid node through tree edge => wrong direction
            rightDir[1] &= false # same as rightDir[1] = false: x && false = false
            DEBUGC && @debug "traverseContainRoot reaches hybrid node $(node.number) through tree edge $(edge.number), so rightDir $(rightDir[1])"
        end
    elseif !node.leaf
        for e in node.edge
            if !isEqual(edge,e) && e.ismajor # minor edges avoided-> their containroot not updated
                other = getOtherNode(e,node);
                if e.containroot # only considered changed those that were true and not hybrid
                    DEBUGC && @debug "traverseContainRoot changing edge $(e.number) to false, at this moment, rightDir is $(rightDir[1])"
                    e.containroot = false;
                    push!(edges_changed, e);
                end
                traverseContainRoot!(other,e, edges_changed, rightDir);
            end
        end
    end
end


# node: must be hybrid node (can come from PhyloNetworks.searchHybridNode)
# return flag, array of edges changed
#        flag: false if the set of edges to place the root is empty
@doc (@doc traverseContainRoot!) updateContainRoot!
function updateContainRoot!(net::HybridNetwork, node::Node)
    node.hybrid || error("node $(node.number )is not hybrid, cannot update containroot")
    edges_changed!(net, Edge[])
    rightDir = [true] #assume good direction, only changed if found hybrid node through tree edge
    for e in node.edge
        if !e.hybrid
            other = getOtherNode(e,node);
            e.containroot = false;
            push!(edges_changed(net),e);
            traverseContainRoot!(other,e, edges_changed(net),rightDir);
        end
    end
    if !rightDir[1] || all((e->!e.containroot), net.edge)
        return false,edges_changed(net)
    else
        return true,edges_changed(net)
    end
end

# function to identify if the network is one of the pathological cases
# see ipad notes: k = 0 (nonidentifiable), k = 1 (nonidentifiable
# "extreme bad triangle", "very" bad triangle I,II) k = 2 (bad diamond
# I,II) also checks if hybrid node has leaf child, in which case,
# major edge is non identifiable
# input: hybrid node around which to check (can come from PhyloNetworks.searchHybridNode)
#        updates gammaz with whatever
# edge lengths are originally in the network
#        allow = true, returns true always, used when reading topology
# returns hasVeryBadTriangle(net), array of edges changed (istIdentifiable, except hybrid edges)
#         false if the network has extremely/very bad triangles
# warning: needs to have updateInCycle already done as it needs inCycle, and k
# check: assume any tree node that has hybrid Edge has only
# one tree edge in cycle (true?)
# updates numBad(net) attribute when found a bad diamond I
function updateGammaz!(net::HybridNetwork, node::Node, allow::Bool)
    node.hybrid || error("node $(node.number) is not hybrid, cannot updategammaz")
    k(node) != -1 || error("update in cycle should have been run before: k(node) not -1")
    isExtBadTriangle!(node, false)
    isVeryBadTriangle!(node, false)
    isBadTriangle!(node, false)
    isBadDiamondI!(node, false)
    isBadDiamondII!(node, false)
    hasVeryBadTriangle!(net, false)
    edges_changed!(net, Edge[])
    edge_maj, edge_min, tree_edge2 = hybridEdges(node);
    other_maj = getOtherNode(edge_maj,node);
    other_min = getOtherNode(edge_min,node);
    k(node) > 2 || error("cycle with only $(k(node)) nodes: parallel edges") # return false, []
    if(k(node) == 4) # could be bad diamond I,II
#        net.numtaxa >= 5 || return false, [] #checked inside optTopRuns now
        edgebla,edge_min2,tree_edge3 = hybridEdges(other_min);
        edgebla,edge_maj2,tree_edge1 = hybridEdges(other_maj);
        other_min2 = getOtherNode(edge_min2,other_min);
        isLeaf1 = getOtherNode(tree_edge1,other_maj);
        isLeaf2 = getOtherNode(tree_edge2,node);
        isLeaf3 = getOtherNode(tree_edge3,other_min);
        tree_edge4 = nothing;
        for e in other_min2.edge
            if(isa(tree_edge4,Nothing) && inCycle(e) == -1 && !e.hybrid)
                tree_edge4 = e;
            end
        end
        if(isEqual(other_min2,getOtherNode(edge_maj2,other_maj)) && isLeaf1.leaf && isLeaf2.leaf && isLeaf3.leaf) # bad diamond I
            @debug "bad diamond I found"
            numBad!(net, numBad(net) + 1)
            isBadDiamondI!(node, true);
            gammaz!(other_min, edge_min.gamma*edge_min2.z);
            gammaz!(other_maj, edge_maj.gamma*edge_maj2.z);
            istIdentifiable!(edge_min2, false)
            istIdentifiable!(edge_maj2, false)
            istIdentifiable!(edge_maj, false)
            istIdentifiable!(edge_min, false)
            push!(edges_changed(net),edge_min2);
            push!(edges_changed(net),edge_min);
            push!(edges_changed(net),edge_maj2);
            push!(edges_changed(net),edge_maj);
        elseif(isEqual(other_min2,getOtherNode(edge_maj2,other_maj)) && isLeaf1.leaf && !isLeaf2.leaf && isLeaf3.leaf && getOtherNode(tree_edge4,other_min2).leaf) # bad diamond II
            @debug "bad diamond II found"
            isBadDiamondII!(node, true);
            setLength!(edge_maj,edge_maj.length+tree_edge2.length)
            setLength!(tree_edge2,0.0)
            push!(edges_changed(net),tree_edge2)
            istIdentifiable!(tree_edge2, false)
            istIdentifiable!(edge_maj, true)
            istIdentifiable!(edge_min, true)
        end
    elseif(k(node) == 3) # could be extreme/very bad triangle or just bad triangle
        if net.numtaxa <= 5
            @debug "extremely or very bad triangle found"
            isVeryBadTriangle!(node, true)
            hasVeryBadTriangle!(net, true)
        elseif net.numtaxa >= 6
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            edgebla,edgebla,tree_edge3 = hybridEdges(other_maj);
            isLeaf1 = getOtherNode(tree_edge1,other_min);
            isLeaf2 = getOtherNode(tree_edge2,node);
            isLeaf3 = getOtherNode(tree_edge3,other_maj);
            if isLeaf1.leaf || isLeaf2.leaf || isLeaf3.leaf
                nl = count([l.leaf for l in [isLeaf1,isLeaf2,isLeaf3]])
                if nl >= 2
                    @debug "warning: extremely bad triangle found"
                    isExtBadTriangle!(node, true);
                    hasVeryBadTriangle!(net, true)
                elseif nl == 1
                    @debug "warning: bad triangle I or II found"
                    isVeryBadTriangle!(node, true);
                    hasVeryBadTriangle!(net, true)
                end
            else
                isBadTriangle!(node, true)
                setLength!(edge_maj,edge_maj.length+tree_edge2.length)
                setLength!(tree_edge2,0.0)
                istIdentifiable!(tree_edge2, false)
                push!(edges_changed(net), tree_edge2);
            end
        end
    end #ends the search for bad things
    if(k(node) > 3 && !isBadDiamondI(node) && !isBadDiamondII(node))
        #println("si entra el ultimo if de k>3 y no bad diamondI,II")
        edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
        if(!istIdentifiable(tree_edge_incycle))
            istIdentifiable!(tree_edge_incycle, true)
            push!(edges_changed(net),tree_edge_incycle);
        end
        istIdentifiable!(edge_maj, isEdgeIdentifiable(edge_maj))
        istIdentifiable!(edge_min, isEdgeIdentifiable(edge_min))
    end
    checkIsBadTriangle(node) == hasVeryBadTriangle(net) || error("node $(node.number) is very bad triangle but hasVeryBadTriangle(net) is $(hasVeryBadTriangle(net))")
    if(allow)
        return true, edges_changed(net)
    else
        return !hasVeryBadTriangle(net), edges_changed(net)
    end
end

updateGammaz!(net::HybridNetwork, node::Node) = updateGammaz!(net, node, false)

#function to check if edge should be identifiable
#it is not only if followed by leaf, or if a newly converted tree edge
function isEdgeIdentifiable(edge::Edge)
    if(edge.hybrid)
        node = edge.node[edge.ischild1 ? 1 : 2]
        #println("is edge $(edge.number) identifiable, node $(node.number)")
        node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
        major,minor,tree = hybridEdges(node)
        #println("major $(major.number), minor $(minor.number), tree $(tree.number)")
        if(getOtherNode(tree,node).leaf)
            return false
        else
            return true
        end
    else
        if(reduce(&,[!edge.node[1].leaf,!edge.node[2].leaf]))
            if(!edge.node[1].hybrid && !edge.node[2].hybrid && !fromBadDiamondI(edge))
                return true
            elseif(edge.node[1].hybrid || edge.node[2].hybrid)
                ind = edge.node[1].hybrid ? 1 : 2
                if(!isBadDiamondII(edge.node[ind]) && !isBadTriangle(edge.node[ind]))
                    return true
                else
                    return false
                end
            end
        else
            return false
        end
    end
end



# function to update the net.partition attribute along a cycle formed
# by nodesChanged vector (obtained from updateInCycle)
# warning: needs updateInCycle for all hybrids before running this
function updatePartition!(net::HybridNetwork, nodesChanged::Vector{Node})
    if net.numhybrids == 0
        net.partition = Partition[]
    end
    for n in nodesChanged
        if(length(n.edge) == 3) #because we are allowing the root to have only two edges when read from parenthetical format
            edge = nothing
            for e in n.edge
                if(inCycle(e) == -1)
                    edge = e
                end
            end
            !isa(edge,Nothing) || error("one edge in n.edge for node $(n.number) should not be in cycle")
            descendants = [edge]
            cycleNum = [inCycle(nodesChanged[1])]
            getDescendants!(getOtherNode(edge,n),edge,descendants,cycleNum)
            !isempty(descendants) || error("descendants is empty for node $(n.number)")
            @debug "for node $(n.number), descendants are $([e.number for e in descendants]), and cycleNum is $(cycleNum)"
            partition = Partition(cycleNum,descendants)
            if(!isPartitionInNet(net,partition)) #need to check not already added by other hybrid nodes
                push!(net.partition, partition)
            end
        end
    end
end

function choosePartition(net::HybridNetwork)
    all((n->(length(n.edges) == 1)), net.partition) && return 0 #cannot put any hyb
    all((n->(length(n.edges) == 3)), net.partition) && return 0 #can only put very bad triangles
    partition = Int[] #good partitions
    for i in 1:length(net.partition)
        if(length(net.partition[i].edges) > 3)
            push!(partition,i)
        end
    end
    isempty(partition) && return 0
    length(partition) == 1 && return partition[1]
    index1 = round(Integer,rand()*size(partition,1));
    while(index1 == 0 || index1 > length(partition))
        index1 = round(Integer,rand()*size(partition,1));
    end
    @debug "chosen partition $([n.number for n in net.partition[partition[index1]].edges])"
    return partition[index1]
end


"""
updateUninformativeQuartets(obsCF::Array{Float64}, tol::Float64)
Returns value to place in quartet.sampled
"""
function updateUninformativeQuartets(quartet::Quartet, atol::Float64)
    allcomp = [(abs(a-b)<=atol) for a in quartet.obsCF', b in quartet.obsCF]
    if false in allcomp
        return(true)
    else
        return(false)
    end
end

"""
updateUninformativeQuartets!(quartets::Vector{Quartet}, tol::Float64)
Checks for quartets classified as 'uninformative' given a tolerance to define CF equality 
    if all CFs in quartet.obsCF are equal within 'tol' tolerance, 
    quartet.sampled is set to 'false'
Output: Updates quartet.sampled in-place for all quartets 
"""
function updateUninformativeQuartets!(quartets::Vector{Quartet}, atol::Float64)
    i = Threads.Atomic{Int}(0);
    Threads.@threads for q in quartets
        q.sampled = updateUninformativeQuartets(q, atol)
        if !(q.sampled)
            #println("bad")
            q.uninformative = true
            Threads.atomic_add!(i, 1)
        end
    end
    return(i[])
end
updateUninformativeQuartets!(d::DataCF, atol::Float64) = updateUninformativeQuartets!(d.quartet, atol)


"""
updateSampledQuartetsAll!(quartets::Vector{Quartet}, toset::Bool)
Updates in-place all Quartet.sampled values where Quartet.uninformative=false 
    to the value provided in toset::Bool
"""
function updateSampledQuartetsAll!(quartets::Vector{Quartet}, toset::Bool)
    Threads.@threads for q in quartets
        if !q.uninformative
            q.sampled=toset
        end
    end
end
updateSampledQuartetsAll!(quartets::Vector{Quartet}) = updateSampledQuartetsAll!(quartets, true)
updateSampledQuartetsAll!(d::DataCF, toset::Bool) = updateSampledQuartetsAll!(d.quartet, true)
updateSampledQuartetsAll!(d::DataCF) = updateSampledQuartetsAll!(d.quartet, true)

"""
    updateSubsetQuartets!(quartets::Vector{Quartet}, prop::Float64, update_uninformative=false)
    updateSubsetQuartets!(d::DataCF, prop::Float64, update_uninformative=false)

Randomly sample a proportion `prop` of quartets, and update in-place
`q.sampled` to true for each quartet `q`.
If `update_uninformative` is false, then we instead sample
a proportion `prop` of *informative* quartets, that is, only sampling from
quartets `q` such that `q.uninformative` is false).

Output: None, `quartets` (or `d.quartet`) is modified in-place
"""
function updateSubsetQuartets!(
    quartets::Vector{Quartet},
    prop::Float64,
    update_uninformative::Bool=false
)
    #reset Quartet.samples for all informative quartets
    updateSampledQuartetsAll!(quartets, false)

    #get count of informative quartets
    indices = [i for (i, q) in enumerate(quartets) if (update_uninformative || !q.uninformative)]
    numQuartets = size(indices, 1)

    #randomly sample numQuartets*propQuartets quartets
    to_sample = Int(round(numQuartets*prop))
    sampled = prop == 1.0 ? indices : sample(indices, to_sample, replace=false)
    for s in sampled 
        quartets[s].sampled = true
    end
end
updateSubsetQuartets!(d::DataCF, prop::Float64,update_uninformative::Bool=false) = updateSubsetQuartets!(d.quartet, prop, update_uninformative)
