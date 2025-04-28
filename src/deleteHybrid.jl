# functions to delete hybridization
# originally in functions.jl
# Claudia MArch 2015

# --------------------------------- delete hybridization -------------------------------

# function to identify inCycle (with priority queue)
# based on updateInCycle as it is the exact same code
# only without changing the incycle attribute
# only returning array of edges/nodes affected by the hybrid
# used when attempting to delete
# input: hybrid node around which we want to identify inCycle
# returns tuple: nocycle, array of edges changed, array of nodes changed
# check: is this traversal much faster than a simple loop over
#        all edges/nodes and check if incycle==hybrid.number?
function identifyInCycle(net::Network,node::Node)
    node.hybrid || error("node $(node.number) is not hybrid, cannot identifyInCycle")
    start = node;
    hybedge = getparentedgeminor(node)
    lastnode = getOtherNode(hybedge,node)
    dist = 0;
    queue = PriorityQueue();
    path = Node[];
    edges_changed!(net, Edge[])
    nodes_changed!(net, Node[])
    push!(edges_changed(net),hybedge);
    push!(nodes_changed(net),node);
    found = false;
    visited!(net, [false for i = 1:size(net.node,1)]);
    enqueue!(queue,node,dist);
    while(!found)
        if(isempty(queue))
            return true, edges_changed(net), nodes_changed(net)
        else
            curr = dequeue!(queue);
            if isEqual(curr,lastnode)
                found = true;
                push!(path,curr);
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
        end
    end # end while
    curr = pop!(path);
    while(!isEqual(curr, start))
        if(inCycle(curr) == start.number)
            push!(nodes_changed(net), curr);
            edge = getconnectingedge(curr, curr.prev);
            if(inCycle(edge) == start.number)
                push!(edges_changed(net), edge);
            end
            curr = curr.prev;
        end
    end
    return false, edges_changed(net), nodes_changed(net)
end



# aux function to traverse the network
# similar to traverseContainRoot but only
# identifying the edges that would be changed by a given
# hybridization
# warning: it does not go accross hybrid node/edge,
#          nor tree node with minor hybrid edge
function traverseIdentifyRoot(node::Node, edge::Edge, edges_changed::Array{Edge,1})
    if(!node.leaf && !node.hybrid)
        for e in node.edge
            if !isEqual(edge,e) && e.ismajor && !e.hybrid
                other = getOtherNode(e,node);
                push!(edges_changed, e);
                if(!other.hybrid)
                    #if(!hasHybEdge(other))
                        traverseIdentifyRoot(other,e, edges_changed);
                    #else
                    #    if hybridEdges(other)[1].ismajor
                    #        traverseIdentifyRoot(other,e, edges_changed);
                    #    end
                    #end
                end
            end
        end
    end
end


# function to identify containroot
# depending on a hybrid node on the network
# input: hybrid node (can come from searchHybridNode)
# return array of edges affected by the hybrid node
function identifyContainRoot(net::HybridNetwork, node::Node)
    node.hybrid || error("node $(node.number) is not hybrid, cannot identify containroot")
    edges_changed!(net, Edge[])
    for e in node.edge
        if(!e.hybrid)
            other = getOtherNode(e,node);
            push!(edges_changed(net),e);
            traverseIdentifyRoot(other,e, edges_changed(net));
        end
    end
    return edges_changed(net)
end

# function to undo the effect of a hybridization
# and then delete it
# input: network, hybrid node, random flag
#        random = true, deletes one hybrid egde at random
#        (minor with prob 1-gamma, major with prob gamma)
#        random = false, deletes the minor edge always
# warning: it uses the gamma of the hybrid edges even if
#          it is not identifiable like in bad diamond I (assumes undone by now)
# blacklist = true: add the edge as a bad choice to put a hybridization (not fully tested)
function deleteHybridizationUpdate!(net::HybridNetwork, hybrid::Node, random::Bool, blacklist::Bool)
    hybrid.hybrid || error("node $(hybrid.number) is not hybrid, so we cannot delete hybridization event around it")
    @debug "MOVE: delete hybridization on hybrid node $(hybrid.number)"
    nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
    !nocycle || error("the hybrid node $(hybrid.number) does not create a cycle")
    edgesRoot = identifyContainRoot(net,hybrid);
    edges = hybridEdges(hybrid);
    undoGammaz!(hybrid,net);
    othermaj = getOtherNode(edges[1],hybrid)
    edgesmaj = hybridEdges(othermaj)
    @debug "edgesmaj[3] $(edgesmaj[3].number) is the one to check if containroot=false already: $(edgesmaj[3].containroot)"
    if edgesmaj[3].containroot #if containroot=true, then we need to undo
        push!(edgesRoot, edges[1]) ## add hybrid edges to edgesRoot to undo containroot
        push!(edgesRoot, edges[2])
        undoContainRoot!(edgesRoot);
    end
    limit = edges[1].gamma
    @debug (limit < 0.5 ? "strange major hybrid edge $(edges[1].number) with γ $limit < 0.5" :
            (limit == 1.0 ? "strange major hybrid edge $(edges[1].number) with γ = $limit" : ""))
    if(random)
        minor = rand() < limit ? false : true
    else
        minor = true;
    end
    deleteHybrid!(hybrid,net,minor, blacklist)
    undoInCycle!(edgesInCycle, nodesInCycle); #moved after deleteHybrid to mantain who is incycle when deleteEdge and look for partition
    undoPartition!(net,hybrid, edgesInCycle)
end

deleteHybridizationUpdate!(net::HybridNetwork, hybrid::Node) = deleteHybridizationUpdate!(net, hybrid, true, false)

# function to delete a hybridization event
# input: hybrid node and network
#        minor: true (deletes minor edge), false (deletes major)
# warning: it is meant after undoing the effect of the
#          hybridization in deleteHybridizationUpdate!
#          by itself, it leaves things as is
# branch lengths of -1.0 are interpreted as missing.
function deleteHybrid!(node::Node,net::HybridNetwork,minor::Bool, blacklist::Bool)
    node.hybrid || error("node $(node.number) has to be hybrid for deleteHybrid")
    if(minor)
        hybedge1,hybedge2,treeedge1 = hybridEdges(node);
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        other3 =  getOtherNode(treeedge1,node);
        if(hybedge1.number > treeedge1.number)
            setLength!(treeedge1, addBL(treeedge1.length, hybedge1.length));
            removeNode!(node,treeedge1);
            setNode!(treeedge1,other1);
            setEdge!(other1,treeedge1);
            removeEdge!(other1, hybedge1);
            PhyloNetworks.deleteEdge!(net,hybedge1);
            #treeedge1.containroot = (!treeedge1.containroot || !hybedge1.containroot) ? false : true #causes problems if hybrid.CR=false
            if(blacklist)
                println("put in blacklist edge $(treeedge1.number)")
                push!(blacklist(net), treeedge1.number)
            end
        else
            makeEdgeTree!(hybedge1,node)
            hasHybEdge!(other1, false);
            setLength!(hybedge1, addBL(hybedge1.length, treeedge1.length));
            removeNode!(node,hybedge1);
            setNode!(hybedge1,other3);
            setEdge!(other3,hybedge1);
            removeEdge!(other3,treeedge1);
            PhyloNetworks.deleteEdge!(net,treeedge1);
            hybedge1.containroot = (!treeedge1.containroot || !hybedge1.containroot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(hybedge1.number)")
                push!(blacklist(net), hybedge1.number)
            end
        end
        hybindex = findfirst([e.hybrid for e in other2.edge]);
        isnothing(hybindex) && error("didn't find hybrid edge in other2")
        if(hybindex == 1)
            treeedge1 = other2.edge[2];
            treeedge2 = other2.edge[3];
        elseif(hybindex == 2)
            treeedge1 = other2.edge[1];
            treeedge2 = other2.edge[3];
        elseif(hybindex == 3)
            treeedge1 = other2.edge[1];
            treeedge2 = other2.edge[2];
        else
            error("strange node has more than three edges")
        end
        treenode1 = getOtherNode(treeedge1,other2);
        treenode2 = getOtherNode(treeedge2,other2);
        if(abs(treeedge1.number) > abs(treeedge2.number))
            setLength!(treeedge2, addBL(treeedge2.length, treeedge1.length));
            removeNode!(other2,treeedge2);
            setNode!(treeedge2,treenode1);
            setEdge!(treenode1,treeedge2);
            removeEdge!(treenode1,treeedge1);
            PhyloNetworks.deleteEdge!(net,treeedge1);
            treeedge2.containroot = (!treeedge1.containroot || !treeedge2.containroot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(treeedge2.number)")
                push!(blacklist(net), treeedge2.number)
            end
        else
            setLength!(treeedge1, addBL(treeedge2.length, treeedge1.length));
            removeNode!(other2,treeedge1);
            setNode!(treeedge1,treenode2);
            setEdge!(treenode2,treeedge1);
            removeEdge!(treenode2,treeedge2);
            PhyloNetworks.deleteEdge!(net,treeedge2);
            treeedge1.containroot = (!treeedge1.containroot || !treeedge2.containroot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(treeedge1.number)")
                push!(blacklist(net), treeedge1.number)
            end
        end
        #removeHybrid!(net,node);
        PhyloNetworks.deleteNode!(net,node);
        PhyloNetworks.deleteNode!(net,other2);
        PhyloNetworks.deleteEdge!(net,hybedge2);
    else
        hybedge1,hybedge2,treeedge1 = hybridEdges(node);
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        setLength!(treeedge1, addBL(treeedge1.length, hybedge2.length))
        removeEdge!(other2,hybedge2)
        removeNode!(node,treeedge1)
        setEdge!(other2,treeedge1)
        setNode!(treeedge1,other2)
        #removeHybrid!(net,node)
        PhyloNetworks.deleteNode!(net,node)
        PhyloNetworks.deleteEdge!(net,hybedge1)
        PhyloNetworks.deleteEdge!(net,hybedge2)
        removeEdge!(other1,hybedge1)
        size(other1.edge,1) == 2 || error("strange node $(other1.number) had 4 edges")
        if(abs(other1.edge[1].number) < abs(other1.edge[2].number))
            edge = other1.edge[1]
            otheredge = other1.edge[2]
        else
            edge = other1.edge[2]
            otheredge = other1.edge[1]
        end
        setLength!(other1.edge[1], addBL(other1.edge[1].length, other1.edge[2].length))
        other3 =  getOtherNode(otheredge,other1);
        removeNode!(other1,edge)
        removeEdge!(other3,otheredge)
        setEdge!(other3,edge)
        setNode!(edge,other3)
        PhyloNetworks.deleteNode!(net,other1)
        PhyloNetworks.deleteEdge!(net,otheredge)
    end
end

deleteHybrid!(node::Node,net::HybridNetwork,minor::Bool) = deleteHybrid!(node,net,minor, false)


# function to update net.partition after deleting a hybrid node
# needs a list of the edges in cycle
function undoPartition!(net::HybridNetwork, hybrid::Node, edgesInCycle::Vector{Edge})
    hybrid.hybrid || error("node $(hybrid.number) is not hybrid, and we need hybrid node inside deleteHybUpdate for undoPartition")
    if net.numhybrids == 0
        net.partition = Partition[]
    else
        cycles = Int[]
        edges = Edge[]
        N = length(net.partition)
        i = 1
        while(i <= N)
            @debug "hybrid number is $(hybrid.number) and partition is $([e.number for e in net.partition[i].edges]), with cycle $(net.partition[i].cycle)"
            if(in(hybrid.number,net.partition[i].cycle))
                @debug "hybrid number matches with partition.cycle"
                p = splice!(net.partition,i)
                @debug "after splice, p partition has edges $([e.number for e in p.edges]) and cycle $(p.cycle)"
                ind = findfirst(isequal(hybrid.number), p.cycle)
                isnothing(ind) && error("hybrid not found in p.cycle")
                deleteat!(p.cycle,ind) #get rid of that hybrid number
                cycles = vcat(cycles,p.cycle)
                edges = vcat(edges,p.edges)
                @debug "edges is $([e.number for e in edges]) and cycles is $(cycles)"
                N = length(net.partition)
            else
                i += 1
            end
        end
        for e in edgesInCycle
            @debug "edge in cycle is $(e.number)"
            if(isEdgeNumIn(e,net.edge)) #only include edge if still in net
                @debug "edge is in net still"
                push!(edges,e)
            end
        end
        newPartition = Partition(unique(cycles),edges)
        @debug "new partition with cycle $(newPartition.cycle), edges $([e.number for e in newPartition.edges])"
        push!(net.partition,newPartition)
    end
end

