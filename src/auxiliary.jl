# auxiliary functions for all the other methods
# originally in functions.jl
# Claudia February 2015
#####################

function setCHECKNET(b::Bool)
    global CHECKNET
    CHECKNET = b
    CHECKNET && @warn "PhyloNetworks.CHECKNET is true: will slow snaq! down."
    b || @info "PhyloNetworks.CHECKNET set to false"
end

# ----- aux general functions ---------------

#based in coupon's collector: E+sqrt(V)
function coupon(n::Number)
    return n*log(n) + n
end

function binom(n::Number,k::Number)
    n >= k || return 0
    n == 1 && return 1
    k == 0 && return 1
    binom(n-1,k-1) + binom(n-1,k) #recursive call
end



"""
    deleteNode!(net::QuartetNetwork, n::Node)

Delete node `n` from a network, i.e. removes it from
net.node, and from net.hybrid or net.leaf as appropriate.
Update attributes `numNodes`, `numTaxa`, `numHybrids`.
Warning: `net.names` is *not* updated, and this is a feature (not a bug)
for networks of type QuartetNetwork.

Warning: if the root is deleted, the new root is arbitrarily set to the
first node in the list. This is intentional to save time because this function
is used frequently in snaq!, which handles semi-directed (unrooted) networks.
"""

function deleteNode!(net::QuartetNetwork, n::Node)
    index = findfirst(no -> no.number == n.number, net.node)
    # isEqual (from above) checks for more than node number
    index !== nothing || error("Node $(n.number) not in quartet network");
    deleteat!(net.node,index);
    net.numNodes -= 1
    if n.hybrid
       removeHybrid!(net,n)
    end
    if n.leaf
        index = findfirst(no -> no === n, net.leaf)
        index !== nothing || error("node $(n.number) not net.leaf")
        deleteat!(net.leaf,index)
        net.numTaxa -= 1
    end
end

"""
    deleteEdge!(net::QuartetNetwork, e::Edge)

Delete edge `e` from `net.edge` and update `net.numEdges`.
If `part` is true, update the network's partition field.
"""


# function to delete an Edge in net.edge and
# update numEdges from a QuartetNetwork
function deleteEdge!(net::QuartetNetwork, e::Edge)
    index = findfirst(x -> x.number == e.number, net.edge)
    # isEqual (from above) checks for more than edge number
    index !== nothing || error("edge not in quartet network");
    deleteat!(net.edge,index);
    net.numEdges -= 1;
end


"""
    removeHybrid!(net::Network, n::Node)

Delete a hybrid node `n` from `net.hybrid`, and update `net.numHybrid`.
The actual node `n` is not deleted. It is kept in the full list `net.node`.
"""
function removeHybrid!(net::Network, n::Node)
    n.hybrid || error("cannot delete node $(n.number) from net.hybrid because it is not hybrid")
    i = findfirst(x -> x===n, net.hybrid)
    i !== nothing || error("hybrid node $(n.number) not in the network's list of hybrids");
    deleteat!(net.hybrid, i);
    net.numHybrids -= 1;
end


"""
    printEdges(io::IO, net)

Print information on the edges of a `QuartetNetwork` object
`net`: edge number, numbers of nodes attached to it, edge length, whether it's
a hybrid edge, its γ inheritance value, whether it's a major edge,
if it could contain the root (this field is not always updated, though)
and attributes pertaining to level-1 networks used in SNaQ:
in which cycle it is contained (-1 if no cycle), and if the edge length
is identifiable (based on quartet concordance factors).
"""


function printEdges(io::IO, net::QuartetNetwork)
    println(io, "edge parent child  length  hybrid isMajor gamma   containRoot inCycle istIdentitiable")
    for e in net.edge
        @printf(io, "%-4d %-6d %-6d ", e.number, getparent(e).number, getchild(e).number)
        @printf(io, "%-7.3f %-6s %-7s ", e.length, e.hybrid, e.isMajor)
        @printf(io, "%-7.4g %-11s %-7d %-5s\n", e.gamma, e.containRoot, e.inCycle, e.istIdentifiable)
    end
end


# ----------------------------------------------------------------------------------------

# setLength
# warning: allows to change edge length for istIdentifiable=false
#          but issues a warning
# negative=true means it allows negative branch lengths (useful in qnet typeHyb=4)
function setLength!(edge::Edge, new_length::Number, negative::Bool)
    (negative || new_length >= 0) || error("length has to be nonnegative: $(new_length), cannot set to edge $(edge.number)")
    new_length >= -0.4054651081081644 || error("length can be negative, but not too negative (greater than -log(1.5)) or majorCF<0: new length is $(new_length)")
    #println("setting length $(new_length) to edge $(edge.number)")
    if(new_length > 10.0)
        new_length = 10.0;
    end
    edge.length = new_length;
    edge.y = exp(-new_length);
    edge.z = 1.0 - edge.y;
    #edge.istIdentifiable || @warn "set edge length for edge $(edge.number) that is not identifiable"
    return nothing
end

"""
    setLength!(edge, newlength)`

Set the length of `edge`, and set `edge.y` and `edge.z` accordingly.
Warning: specific to SNaQ. Use [`setlengths!`](@ref) or [`setBranchLength!`](@ref)
for more general tools.

- The new length is censored to 10: if the new length is above 10,
  the edge's length will be set to 10. Lengths are interpreted in coalescent
  units, and 10 is close to infinity: near perfect gene tree concordance.
  10 is used as an upper limit to coalescent units that can be reliably estimated.
- The new length is allowed to be negative, but must be greater than -log(1.5),
  to ensure that the major quartet concordance factor (1 - 2/3 exp(-length)) is >= 0.
"""
setLength!(edge::Edge, new_length::Number) = setLength!(edge, new_length, false)



"""
    setGammaBLfromGammaz!(node, network)

Update the γ values of the two sister hybrid edges in a bad diamond I, given the `gammaz` values
of their parent nodes, and update the branch lengths t1 and t2 of their parent edges
(those across from the hybrid nodes), in such a way that t1=t2 and that these branch lengths
and γ values are consistent with the `gammaz` values in the network.

Similar to the first section of [`undoGammaz!`](@ref),
but does not update anything else than γ and t's.
Unlike `undoGammaz!`, no error if non-hybrid `node` or not at bad diamond I.
"""
function setGammaBLfromGammaz!(node::Node, net::HybridNetwork)
    if !node.isBadDiamondI || !node.hybrid
        return nothing
    end
    edge_maj, edge_min, tree_edge2 = hybridEdges(node);
    other_maj = getOtherNode(edge_maj,node);
    other_min = getOtherNode(edge_min,node);
    edgebla,tree_edge_incycle1,tree_edge = hybridEdges(other_min);
    edgebla,tree_edge_incycle2,tree_edge = hybridEdges(other_maj);
    if(approxEq(other_maj.gammaz,0.0) && approxEq(other_min.gammaz,0.0))
        edge_maj.gamma = 1.0 # γ and t could be anything if both gammaz are 0
        edge_min.gamma = 0.0 # will set t's to 0 and minor γ to 0.
        newt = 0.0
    else
        ((approxEq(other_min.gammaz,0.0) || other_min.gammaz >= 0.0) &&
         (approxEq(other_maj.gammaz,0.0) || other_maj.gammaz >= 0.0)    ) ||
            error("bad diamond I in node $(node.number) but missing (or <0) gammaz")
        ztotal = other_maj.gammaz + other_min.gammaz
        edge_maj.gamma = other_maj.gammaz / ztotal
        edge_min.gamma = other_min.gammaz / ztotal
        newt = -log(1-ztotal)
    end
    setLength!(tree_edge_incycle1,newt)
    setLength!(tree_edge_incycle2,newt)
end

# function to find if a given partition is in net.partition
function isPartitionInNet(net::HybridNetwork,desc::Vector{Edge},cycle::Vector{Int})
    for p in net.partition
        if(sort(cycle) == sort(p.cycle))
            if(sort([e.number for e in desc]) == sort([e.number for e in p.edges]))
                return true
            end
        end
    end
    return false
end

# function to check if a partition is already in net.partition
# used in updatePartition
function isPartitionInNet(net::HybridNetwork,partition::Partition)
    if(isempty(net.partition))
        return false
    end
    for p in net.partition
        cycle = isempty(setdiff(p.cycle,partition.cycle)) && isempty(setdiff(partition.cycle,p.cycle))
        edges = isempty(setdiff([n.number for n in p.edges],[n.number for n in partition.edges])) && isempty(setdiff([n.number for n in partition.edges],[n.number for n in p.edges]))
        if(cycle && edges)
            return true
        end
    end
    return false
end

# function to check that everything matches in a network
# in particular, cycles, partitions and containRoot
# fixit: need to add check on identification of bad diamonds, triangles
# and correct computation of gammaz
# light=true: it will not collapse with nodes with 2 edges, will return a flag of true
# returns true if found egde with BL -1.0 (only when light=true, ow error)
# added checkPartition for undirectedOtherNetworks that do not need correct hybrid node number
function checkNet(net::HybridNetwork, light::Bool; checkPartition=true::Bool)
    @debug "checking net"
    net.numHybrids == length(net.hybrid) || error("discrepant number on net.numHybrids (net.numHybrids) and net.hybrid length $(length(net.hybrid))")
    net.numTaxa == length(net.leaf) || error("discrepant number on net.numTaxa (net.numTaxa) and net.leaf length $(length(net.leaf))")
    net.numNodes == length(net.node) || error("discrepant number on net.numNodes (net.numNodes) and net.node length $(length(net.node))")
    net.numEdges == length(net.edge) || error("discrepant number on net.numEdges (net.numEdges) and net.edge length $(length(net.edge))")
    if(isTree(net))
        all(x->x.containRoot,net.edge) || error("net is a tree, but not all edges can contain root")
        all(x->x.isMajor,net.edge) || error("net is a tree, but not all edges are major")
        all(x->!(x.hybrid),net.edge) || error("net is a tree, but not all edges are tree")
        all(x->!(x.hybrid),net.node) || error("net is a tree, but not all nodes are tree")
        all(x->!(x.hasHybEdge),net.node) || error("net is a tree, but not all nodes hasHybEdge=false")
        all(x->(x.gamma == 1.0 ? true : false),net.edge) || error("net is a tree, but not all edges have gamma 1.0")
    end
    for h in net.hybrid
        if(isBadTriangle(h))
            @debug "hybrid $(h.number) is very bad triangle"
            net.hasVeryBadTriangle || error("hybrid node $(h.number) is very bad triangle, but net.hasVeryBadTriangle is $(net.hasVeryBadTriangle)")
            h.isVeryBadTriangle || h.isExtBadTriangle || error("hybrid node $(h.number) is very bad triangle but it does not know it")
        end
        nocycle,edges,nodes = identifyInCycle(net,h)
        for e in edges
            e.inCycle == h.number || error("edge $(e.number) is in cycle of hybrid node $(h.number) but its inCycle attribute is $(e.inCycle)")
            if(e.length == -1.0)
                if(light)
                    return true
                else
                    error("found edge with BL -1.0")
                end
            end
            if(e.hybrid)
                !e.containRoot || error("hybrid edge $(e.number) should not contain root") # fixit: disagree
                o = getOtherNode(e,h)
                o.hasHybEdge || error("found node $(o.number) attached to hybrid edge but hasHybEdge=$(o.hasHybEdge)")
            end
        end
        for n in nodes
            n.inCycle == h.number || error("node $(n.number) is in cycle of hybrid node $(h.number) but its inCycle attribute is $(n.inCycle)")
            e1,e2,e3 = hybridEdges(n)
            i = 0
            for e in [e1,e2,e3]
                if(isa(e,Nothing) && h.k != 2)
                    error("edge found that is Nothing, and hybrid node $(h.number) k is $(h.k). edge as nothing can only happen when k=2")
                elseif(!isa(e,Nothing))
                    if(e.inCycle == -1)
                        i += 1
                        desc = [e]
                        cycleNum = [h.number]
                        getDescendants!(getOtherNode(e,n),e,desc,cycleNum)
                        if(checkPartition && !isPartitionInNet(net,desc,cycleNum))
                            printPartitions(net)
                            error("partition with cycle $(cycleNum) and edges $([e.number for e in desc]) not found in net.partition")
                        end
                    end
                end
            end
            i == 1 || error("strange node $(n.number) incycle $(h.number) but with $(i) edges not in cycle, should be only one")
            edgesRoot = identifyContainRoot(net,h)
            for edge in edgesRoot
                if edge.containRoot
                    @debug begin printEverything(net); "printed everything" end
                    error("edge $(edge.number) should not contain root")
                end
            end
        end
    end
    for n in net.node
        if(n.leaf)
            length(n.edge) == 1 || error("leaf $(n.number) with $(length(n.edge)) edges instead of 1")
        else
            if(light)
                if(length(n.edge) != 3)
                    @debug "warning: node $(n.number) with $(length(n.edge)) edges instead of 3"
                    return true
                end
            else
                length(n.edge) == 3 || error("node $(n.number) with $(length(n.edge)) edges instead of 3")
            end
        end
    end
    for e in net.edge
        if(e.length == -1.0)
            if(light)
                return true
            else
                error("edge found with BL -1.0")
            end
        end
    end
    @debug "no errors in checking net"
    return false
end

checkNet(net::HybridNetwork) = checkNet(net, false)


# function to print everything for a given net
# this is used a lot inside snaq to debug, so need to use level1 attributes
# and not change the network: with writeTopologyLevel1
function printEverything(net::HybridNetwork)
    printEdges(net)
    printNodes(net)
    printPartitions(net)
    println("$(writeTopologyLevel1(net))")
end

# function to check if a node is very or ext bad triangle
function isBadTriangle(node::Node)
    node.hybrid || error("cannot check if node $(node.number) is very bad triangle because it is not hybrid")
    if(node.k == 3)
        edgemaj, edgemin, treeedge = hybridEdges(node)
        othermaj = getOtherNode(edgemaj,node)
        othermin = getOtherNode(edgemin,node)
        treenode = getOtherNode(treeedge,node)
        edges1 = hybridEdges(othermaj)
        o1 = getOtherNode(edges1[3],othermaj)
        edges2 = hybridEdges(othermin)
        o2 = getOtherNode(edges2[3],othermin)
        leaves = sum([n.leaf ? 1 : 0 for n in [treenode,o1,o2]])
        if(leaves == 1 || leaves == 2)
            return true
        else
            return false
        end
    else
        return false
    end
end


# function to switch a hybrid node in a network to another node in the cycle
function switchHybridNode!(net::HybridNetwork, hybrid::Node, newHybrid::Node)
    hybrid.hybrid || error("node $(hybrid.number) has to be hybrid to switch to a different hybrid")
    newHybrid.inCycle == hybrid.number || error("new hybrid needs to be in the cycle of old hybrid: $(hybrid.number)")
    !newHybrid.hybrid || error("strange hybrid node $(newHybrid.number) in cycle of another hybrid $(hybrid.number)")
    newHybrid.hybrid = true
    newHybrid.hasHybEdge = true
    newHybrid.name = hybrid.name
    pushHybrid!(net,newHybrid)
    makeNodeTree!(net,hybrid)
end

"""
    sorttaxa!(DataFrame, columns)

Reorder the 4 taxa and reorders the observed concordance factors accordingly, on each row of
the data frame. If `columns` is ommitted, taxon names are assumed to be in columns 1-4 and
CFs are assumed to be in columns 5-6 with quartets in this order: `12_34`, `13_24`, `14_23`.
Does **not** reorder credibility interval values, if present.

    sorttaxa!(DataCF)
    sorttaxa!(Quartet, permutation_tax, permutation_cf)

Reorder the 4 taxa in each element of the DataCF `quartet`. For a given Quartet,
reorder the 4 taxa in its fields `taxon` and `qnet.quartetTaxon` (if non-empty)
and reorder the 3 concordance values accordingly, in `obsCF` and `qnet.expCF`.

`permutation_tax` and `permutation_cf` should be vectors of short integers (Int8) of length 4 and 3
respectively, whose memory allocation gets reused. Their length is *not checked*.

`qnet.names` is unchanged: the order of taxon names here relates to the order of nodes in the network
(???)
"""
function sorttaxa!(dat::DataCF)
    ptax = Array{Int8}(undef, 4) # to hold the sort permutations
    pCF  = Array{Int8}(undef, 3)
    for q in dat.quartet
        sorttaxa!(q, ptax, pCF)
    end
end

function sorttaxa!(df::DataFrame, co=Int[]::Vector{Int})
    if length(co)==0
        co = collect(1:7)
    end
    length(co) > 6 || error("column vector must be of length 7 or more")
    ptax = Array{Int8}(undef, 4)
    pCF  = Array{Int8}(undef, 3)
    taxnam = Array{eltype(df[!,co[1]])}(undef, 4)
    for i in 1:size(df,1)
        for j=1:4 taxnam[j] = df[i,co[j]]; end
        sortperm!(ptax, taxnam)
        sorttaxaCFperm!(pCF, ptax) # update permutation pCF according to taxon permutation
        df[i,co[1]], df[i,co[2]], df[i,co[3]], df[i,co[4]] = taxnam[ptax[1]], taxnam[ptax[2]], taxnam[ptax[3]], taxnam[ptax[4]]
        df[i,co[5]], df[i,co[6]], df[i,co[7]] = df[i,co[pCF[1]+4]], df[i,co[pCF[2]+4]], df[i,co[pCF[3]+4]]
    end
    return df
end

function sorttaxa!(qua::Quartet, ptax::Vector{Int8}, pCF::Vector{Int8})
    qt = qua.taxon
    if length(qt)==4
        sortperm!(ptax, qt)
        sorttaxaCFperm!(pCF, ptax) # update permutation pCF accordingly
        qt[1], qt[2], qt[3], qt[4] = qt[ptax[1]], qt[ptax[2]], qt[ptax[3]], qt[ptax[4]]
        qua.obsCF[1], qua.obsCF[2], qua.obsCF[3] = qua.obsCF[pCF[1]], qua.obsCF[pCF[2]], qua.obsCF[pCF[3]]
        # do *NOT* modify qua.qnet.quartetTaxon: it points to the same array as qua.taxon
        eCF = qua.qnet.expCF
        if length(eCF)==3
            eCF[1], eCF[2], eCF[3] = eCF[pCF[1]], eCF[pCF[2]], eCF[pCF[3]]
        end
    elseif length(qt)!=0
        error("Quartet with $(length(qt)) taxa")
    end
    return qua
end

# find permutation pCF of the 3 CF values: 12_34, 13_24, 14_23. 3!=6 possible permutations
# ptax = one of 4!=24 possible permutations on the 4 taxon names
# kernel: pCF = identity if ptax = 1234, 2143, 3412 or 4321
# very long code, but to minimize equality checks at run time
function sorttaxaCFperm!(pcf::Vector{Int8}, ptax::Vector{Int8})
    if ptax[1]==1
        if     ptax[2]==2
            pcf[1]=1
            if  ptax[3]==3 # ptax = 1,2,3,4
                pcf[2]=2; pcf[3]=3
            else           # ptax = 1,2,4,3
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==3
            pcf[1]=2
            if  ptax[3]==2 # ptax = 1,3,2,4
                pcf[2]=1; pcf[3]=3
            else           # ptax = 1,3,4,2
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==4
            pcf[1]=3
            if  ptax[3]==2 # ptax = 1,4,2,3
                pcf[2]=1; pcf[3]=2
            else           # ptax = 1,4,3,2
                pcf[2]=2; pcf[3]=1
            end
        end
    elseif ptax[1]==2
        if     ptax[2]==1
            pcf[1]=1
            if  ptax[3]==4 # ptax = 2,1,4,3
                pcf[2]=2; pcf[3]=3
            else           # ptax = 2,1,3,4
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==4
            pcf[1]=2
            if  ptax[3]==1 # ptax = 2,4,1,3
                pcf[2]=1; pcf[3]=3
            else           # ptax = 2,4,3,1
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==3
            pcf[1]=3
            if  ptax[3]==1 # ptax = 2,3,1,4
                pcf[2]=1; pcf[3]=2
            else           # ptax = 2,3,4,1
                pcf[2]=2; pcf[3]=1
            end
        end
    elseif ptax[1]==3
        if     ptax[2]==4
            pcf[1]=1
            if  ptax[3]==1 # ptax = 3,4,1,2
                pcf[2]=2; pcf[3]=3
            else           # ptax = 3,4,2,1
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==1
            pcf[1]=2
            if  ptax[3]==4 # ptax = 3,1,4,2
                pcf[2]=1; pcf[3]=3
            else           # ptax = 3,1,2,4
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==2
            pcf[1]=3
            if  ptax[3]==4 # ptax = 3,2,4,1
                pcf[2]=1; pcf[3]=2
            else           # ptax = 3,2,1,4
                pcf[2]=2; pcf[3]=1
            end
        end
    else # ptax[1]==4
        if     ptax[2]==3
            pcf[1]=1
            if  ptax[3]==2 # ptax = 4,3,2,1
                pcf[2]=2; pcf[3]=3
            else           # ptax = 4,3,1,2
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==2
            pcf[1]=2
            if  ptax[3]==3 # ptax = 4,2,3,1
                pcf[2]=1; pcf[3]=3
            else           # ptax = 4,2,1,3
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==1
            pcf[1]=3
            if  ptax[3]==3 # ptax = 4,1,3,2
                pcf[2]=1; pcf[3]=2
            else           # ptax = 4,1,2,3
                pcf[2]=2; pcf[3]=1
            end
        end
    end
end

# function to check in an edge is in an array by comparing
# edge numbers (could use isEqual for adding comparisons of gammaz and inCycle)
# needed for updateHasEdge
function isEdgeNumIn(edge::Edge,array::Array{Edge,1})
    enum = edge.number
    return any(e -> e.number == enum, array)
end

# function to check in a leaf is in an array by comparing
# the numbers (uses isEqual)
# needed for updateHasEdge
function isNodeNumIn(node::Node,array::Array{Node,1})
    return all((e->!isEqual(node,e)), array) ? false : true
end



#------------------------------------
function citation()
    bibfile = joinpath(@__DIR__, "..", "CITATION.bib")
    out = readlines(bibfile)
    println("Bibliography in bibtex format also in CITATION.bib")
    println(join(out,'\n'))
end
