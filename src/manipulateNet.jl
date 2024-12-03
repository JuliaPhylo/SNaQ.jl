"""
    undirectedOtherNetworks(net::HybridNetwork)

Return a vector of HybridNetwork objects, obtained by switching the placement
of each hybrid node to other nodes inside its cycle. This amounts to changing
the direction of a gene flow event (recursively to move around the whole cycle
of each reticulation).

Optional argument: `outgroup`, as a String. If an outgroup is specified,
then networks conflicting with the placement of the root are avoided.

Assumptions: `net` is assumed to be of level 1, that is, each blob has a
single cycle with a single reticulation.
All level-1 fields of `net` are assumed up-to-date.
# Example
```julia
julia> net = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
julia> vnet = undirectedOtherNetworks(net)
```
"""
function undirectedOtherNetworks(net0::HybridNetwork; outgroup="none"::AbstractString, insideSnaq=false::Bool)
# extra optional argument: "insideSnaq". When true, all level1-attributes are assumed up-to-date
# So far, undirectedOtherNetworks is called inside optTopRuns only
# Potential bug: if new node is -1, then inCycle will become meaningless: changed in readSubTree here
# WARNING: does not update partition, because only thing to change is hybrid node number
    if !insideSnaq
        net0 = readnewick_level1(writenewick_level1(net0))
    end
    otherNet = HybridNetwork[]
    for i in 1:net0.numhybrids #need to do for by number, not node
        net = deepcopy(net0) # to avoid redoing attributes after each cycle is finished
        ## undo attributes at current hybrid node:
        hybrid = net.hybrid[i]
        nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
        @debug "nodesInCycle are: $([n.number for n in nodesInCycle])"
        !nocycle || error("the hybrid node $(hybrid.number) does not create a cycle")
        edgesRoot = identifyContainRoot(net,hybrid);
        edges = hybridEdges(hybrid);
        undoGammaz!(hybrid,net);
        othermaj = getOtherNode(edges[1],hybrid)
        edgesmaj = hybridEdges(othermaj)
        if edgesmaj[3].containroot #if containroot=true, then we need to undo
            undoContainRoot!(edgesRoot);
        end
        ## changes to new hybrid node:
        for newn in nodesInCycle
            if newn.number != hybrid.number # nodesInCycle contains the hybrid too
                newnet = deepcopy(net)
                newnocycle, newedgesInCycle, newnodesInCycle = identifyInCycle(newnet,newnet.hybrid[i]);
                !newnocycle || error("the hybrid node $(newnet.hybrid[i].number) does not create a cycle")
                ind = getIndexNode(newn.number,newnet) # find the newn node in the new network
                @debug "moving hybrid to node $(newnet.node[ind].number)"
                hybridatnode!(newnet, newnet.hybrid[i], newnet.node[ind])
                @debug begin printedges(newnet); "printed edges" end
                @debug begin printnodes(newnet); "printed nodes" end
                undoInCycle!(newedgesInCycle, newnodesInCycle);
                @debug begin printedges(newnet); "printed edges" end
                @debug begin printnodes(newnet); "printed nodes" end
                ##undoPartition!(net,hybrid, edgesInCycle)
                success, hybrid0, flag, nocycle, flag2, flag3 = updateAllNewHybrid!(newnet.node[ind], newnet, false,false,false)
                if success
                    @debug "successfully added new network: $(writenewick_level1(newnet))"
                    push!(otherNet,newnet)
                else
                    println("the network obtained by putting the new hybrid in node $(newnet.node[ind].number) is not good, inCycle,gammaz,containroot: $([flag,flag2,flag3]), we will skip it")
                end
            end
        end
    end
    # check root in good position
    if outgroup == "none"
        for n in otherNet
            !isTree(n) && checkRootPlace!(n, verbose=false)
        end
        return otherNet
    else ## root already in good place
        @debug "we will remove networks contradicting the outgroup in undirectedOtherNetworks"
        whichKeep = ones(Bool,length(otherNet)) # repeats 'true'
        i = 1
        for n in otherNet
            if !isTree(n)
                try
                    checkRootPlace!(n, verbose=true, outgroup=outgroup)
                catch
                    @debug "found one network incompatible with outgroup"
                    @debug "$(writenewick_level1(n))"
                    whichKeep[i] = false
                end
            end
            i = i+1;
        end
        return otherNet[whichKeep]
    end
end

"""
    hybridatnode!(net::HybridNetwork, nodeNumber::Integer)

Change the direction and status of edges in network `net`,
to move the hybrid node in a cycle to the node with number `nodeNumber`.
This node must be in one (and only one) cycle, otherwise an error will be thrown.
Check and update the nodes' field `inCycle`.

Output: `net` after hybrid modification.

Assumption: `net` must be of level 1, that is, each blob has a
single cycle with a single reticulation.

# example

```julia
net = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
using PhyloPlots
plot(net, shownodenumber=true); # to locate nodes and their numbers. D of hybrid origin
hybridatnode!(net, -4)
plot(net, shownodenumber=true); # hybrid direction reversed: now 2B of hybrid origin
```
"""
function hybridatnode!(net::HybridNetwork, nodeNumber::Integer)
    undoInCycle!(net.edge, net.node)
    for n in net.hybrid
        flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
        flag || error("not level1 network, hybrid $(n.number) cycle intersects another cycle")
        !nocycle || error("strange network without cycle for hybrid $(n.number)")
    end
    ind = 0
    try
        ind = getIndexNode(nodeNumber,net)
    catch
        error("cannot set node $(nodeNumber) as hybrid because it is not part of net")
    end
    inCycle(net.node[ind]) != -1 || error("node $(nodeNumber) is not part of any cycle, so we cannot make it hybrid")
    indhyb = 0
    try
        indhyb = getIndexNode(inCycle(net.node[ind]),net)
    catch
        error("cannot find the hybrid node with number $(inCycle(net.node[ind]))")
    end
    hybrid = net.node[indhyb]
    hybridatnode!(net,hybrid,net.node[ind])
    return net
end

"""
    hybridatnode!(net::HybridNetwork, hybrid::Node, newNode::Node)

Move the reticulation from `hybrid` to `newNode`,
which must in the same cycle. `net` is assumed to be of level 1,
but **no checks** are made and fields are supposed up-to-date.

Called by `hybridatnode!(net, nodenumber)`, which is itself
called by [`undirectedOtherNetworks`](@ref).
"""
function hybridatnode!(net::HybridNetwork, hybrid::Node, newNode::Node)
    hybrid.hybrid || error("node $(hybrid.number) should be hybrid, but it is not")
    hybedges = hybridEdges(hybrid)
    makeEdgeTree!(hybedges[1],hybrid)
    makeEdgeTree!(hybedges[2],hybrid)
    inCycle!(hybedges[1], hybrid.number) #just to keep attributes ok
    inCycle!(hybedges[2], hybrid.number)
    switchHybridNode!(net,hybrid,newNode)
    found = false
    for e in newNode.edge
        if inCycle(e) == hybrid.number
            if !found
                found = true
                makeEdgeHybrid!(e,newNode, 0.51, switchHyb=true) #first found, major edge, need to optimize gamma anyway
                ##e.gamma = -1
                ##e.containroot = true ## need attributes like in snaq
            else
                makeEdgeHybrid!(e,newNode, 0.49, switchHyb=true) #second found, minor edge
                ##e.gamma = -1
                ##e.containroot = true
            end
        end
    end
end

# Not used anywhere, but tested
# does not call hybridatnode! but repeats its code: oops! violates DRY principle
# nodeNumber should correspond to the number assigned by readnewick_level1,
# and the node numbers in `net` are irrelevant.
"""
    hybridatnode(net::HybridNetwork, nodeNumber::Integer)

Move the hybrid node in a cycle to make node number `nodeNumber` a hybrid node
Compared to [`hybridatnode!`], this method checks that `net` is of level 1
(required) and does not modify it.
"""
function hybridatnode(net0::HybridNetwork, nodeNumber::Integer)
    net = readnewick_level1(writenewick_level1(net0)) # we need inCycle attributes
    ind = 0
    try
        ind = getIndexNode(nodeNumber,net)
    catch
        error("cannot set node $(nodeNumber) as hybrid because it is not part of net")
    end
    inCycle(net.node[ind]) != -1 || error("node $(nodeNumber) is not part of any cycle, so we cannot make it hybrid")
    indhyb = 0
    try
        indhyb = getIndexNode(inCycle(net.node[ind]),net)
    catch
        error("cannot find the hybrid node with number $(inCycle(net.node[ind]))")
    end
    hybrid = net.node[indhyb]
    hybrid.hybrid || error("node $(hybrid.number) should be hybrid, but it is not")
    hybedges = hybridEdges(hybrid)
    makeEdgeTree!(hybedges[1],hybrid)
    makeEdgeTree!(hybedges[2],hybrid)
    inCycle!(hybedges[1], hybrid.number) #just to keep attributes ok
    inCycle!(hybedges[2], hybrid.number)
    switchHybridNode!(net,hybrid,net.node[ind])
    found = false
    for e in net.node[ind].edge
        if inCycle(e) == hybrid.number
            if !found
                found = true
                makeEdgeHybrid!(e,net.node[ind], 0.51, switchHyb=true) #first found, major edge, need to optimize gamma anyway
                e.gamma = -1
                e.containroot = true
            else
                makeEdgeHybrid!(e,net.node[ind], 0.49, switchHyb=true) #second found, minor edge
                e.gamma = -1
                e.containroot = true
            end
        end
    end
    return net
end