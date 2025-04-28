# test functions for the 5taxon networks
# functions used in tests_5taxon_readTopology.jl
# fixit: it would be better that within each function, it would look
# for the needed edges/nodes starting from the hybrid node
# so that it would be a general function (not bounded to the specific indeces)
# Claudia November 2014

# Case C bad triangle II
function testCaseC(net::HybridNetwork)
    n=searchHybridNode(net);
    node = n[1];
    k(node) != 3 ? error("k diff than 3") : nothing
    isVeryBadTriangle(node) ? nothing : error("does not know it is very bad triangle")
    isExtBadTriangle(node) ? error("thinks it is extremely bad triangle") : nothing
    hasVeryBadTriangle(net) ? nothing : error("net does not know it has very bad triangle")
    numBad(net) == 0 ? nothing : error("numBad(net) should be 0")
    net.numhybrids != 1 ? error("should have 1 hybrid, but net.numhybrids is $(net.numhybrids): $([n.number for n in net.hybrid])") : nothing
end

# Case F bad diamond I
function testCaseF(net::HybridNetwork)
    visited!(net, [istIdentifiable(e) for e in net.edge]);
    net.numhybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numhybrids)")
    node=net.hybrid[1];
    k(node) != 4 ? error("k diff than 4") : nothing
    (inCycle(net.edge[3]) != node.number || inCycle(net.edge[4]) != node.number || inCycle(net.edge[5]) != node.number || inCycle(net.edge[7]) != node.number ) ? error("edges not correctly in cycle") : nothing
    (inCycle(net.node[3])  != node.number || inCycle(net.node[4])  != node.number || inCycle(net.node[6])  != node.number || inCycle(net.node[7])  != node.number) ? error("nodes not correctly in cycle") : nothing
    !isBadDiamondI(node) ? error("does not know it is bad diamond I") : nothing
    isBadDiamondII(node) ? error("thinks it is bad diamond II") : nothing
    net.edge[2].containroot ? error("edge can contain root") : nothing
    (!net.edge[3].hybrid || !net.edge[3].ismajor) ? error("edge 3 is not hybrid or major") : nothing
    (!net.edge[5].hybrid || net.edge[5].ismajor) ? error("edge 5 is not hybrid or is major") : nothing
    gammaz(net.node[4]) != net.edge[3].gamma*net.edge[4].z ? error("node 4 gammaz not correctly calculated") : nothing
    gammaz(net.node[6]) != net.edge[5].gamma*net.edge[7].z ? error("node 6 gammaz not correctly calculated") : nothing
    (istIdentifiable(net.edge[4]) || istIdentifiable(net.edge[7])) ? error("edges identifiable that should not") : nothing
    !istIdentifiable(net.edge[8]) ? error("edge 8 not identifiable") : nothing
    visited(net)[8] = false;
    !all([!id for id in visited(net)]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,5,7,8] ? nothing : error("net.leaf is wrong")
end

# Case G
function testCaseG(net::HybridNetwork)
    visited!(net, [istIdentifiable(e) for e in net.edge]);
    net.numhybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numhybrids)")
    node=net.hybrid[1];
    k(node) != 4 ? error("k diff than 4") : nothing
    (inCycle(net.edge[5]) != node.number || inCycle(net.edge[6]) != node.number || inCycle(net.edge[7]) != node.number || inCycle(net.edge[9]) != node.number ) ? error("edges not correctly in cycle") : nothing
    (inCycle(net.node[5])  != node.number || inCycle(net.node[6])  != node.number || inCycle(net.node[8])  != node.number || inCycle(net.node[9])  != node.number) ? error("nodes not correctly in cycle") : nothing
    (isBadDiamondI(node) || isBadDiamondII(node) ) ? error("thinks it is bad diamond") : nothing
    net.edge[4].containroot ? error("edge can contain root") : nothing
    (!net.edge[5].hybrid || !net.edge[5].ismajor) ? error("edge 5 is not hybrid or major") : nothing
    (!net.edge[7].hybrid || net.edge[7].ismajor) ? error("edge 7 is not hybrid or is major") : nothing
    gammaz(net.node[5]) != -1 ? error("node 5 gammaz should be -1") : nothing
    (!istIdentifiable(net.edge[6]) || !istIdentifiable(net.edge[3]) || !istIdentifiable(net.edge[9])) ? error("edges identifiable as not identifiable") : nothing
    visited(net)[6] = false;
    visited(net)[3] = false;
    visited(net)[9] = false;
    !all([!id for id in visited(net)]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,7,8] ? nothing : error("net.leaf is wrong")
end

# Case H
function testCaseH(net::HybridNetwork)
    visited!(net, [istIdentifiable(e) for e in net.edge]);
    net.numhybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numhybrids)")
    node=net.hybrid[1];
    k(node) != 4 ? error("k diff than 4") : nothing
    (inCycle(net.edge[4]) != node.number || inCycle(net.edge[5]) != node.number || inCycle(net.edge[7]) != node.number || inCycle(net.edge[9]) != node.number ) ? error("edges not correctly in cycle") : nothing
    (inCycle(net.node[4])  != node.number || inCycle(net.node[6])  != node.number || inCycle(net.node[8])  != node.number || inCycle(net.node[10])  != node.number) ? error("nodes not correctly in cycle") : nothing
    (isBadDiamondI(node) || isBadDiamondII(node)) ? error("thinks it is bad diamond") : nothing
    net.edge[8].containroot ? error("8 can contain root") : nothing
    (!net.edge[9].hybrid || !net.edge[9].ismajor) ? error("edge 9 is not hybrid or major") : nothing
    (!net.edge[4].hybrid || net.edge[4].ismajor) ? error("edge 4 is not hybrid or is major") : nothing
    gammaz(node) != -1 ? error("hybrid node gammaz should be -1") : nothing
    (!istIdentifiable(net.edge[3]) || !istIdentifiable(net.edge[5]) || !istIdentifiable(net.edge[7])) ? error("edge9,5,13not identifiable") : nothing
    visited(net)[3] = false;
    visited(net)[5] = false;
    visited(net)[7] = false;
    !all([!id for id in visited(net)]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,5,6] ? nothing : error("net.leaf is wrong")
end


# Case J
function testCaseJ(net::HybridNetwork)
    visited!(net, [istIdentifiable(e) for e in net.edge]);
    net.numhybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numhybrids)")
    node=net.hybrid[1];
    k(node) != 5 ? error("k diff than 5") : nothing
    (inCycle(net.edge[2]) != node.number || inCycle(net.edge[4]) != node.number || inCycle(net.edge[6]) != node.number || inCycle(net.edge[10]) != node.number || inCycle(net.edge[8]) != node.number) ? error("edges not correctly in cycle") : nothing
    (inCycle(net.node[2])  != node.number || inCycle(net.node[4])  != node.number || inCycle(net.node[6])  != node.number || inCycle(net.node[10])  != node.number || inCycle(net.node[9])) != node.number ? error("nodes not correctly in cycle") : nothing
    net.edge[1].containroot ? error("edge can contain root") : nothing
    (!net.edge[2].hybrid || !net.edge[2].ismajor) ? error("edge 2 is not hybrid or major") : nothing
    gammaz(node) != -1 ? error("hybrid node gammaz should be -1") : nothing
    (!istIdentifiable(net.edge[4]) || !istIdentifiable(net.edge[6]) || !istIdentifiable(net.edge[10])) ? error("edge9,5,10not identifiable") : nothing
    visited(net)[4] = false;
    visited(net)[6] = false;
    visited(net)[10] = false;
    !all([!id for id in visited(net)]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,3,4,5,6] ? nothing : error("net.leaf is wrong")
end

# Case D bad triangle I
function testCaseD(net::HybridNetwork)
 visited!(net, [istIdentifiable(e) for e in net.edge]);
    n = searchHybridNode(net);
    node = n[1];
    k(node) != 3 ? error("k diff than 3") : nothing
    isVeryBadTriangle(node) ? nothing : error("does not know it is very bad triangle")
    isExtBadTriangle(node) ? error("thinks it is extremely bad triangle") : nothing
    hasVeryBadTriangle(net) ? nothing : error("net does not know it has very bad triangle")
    numBad(net) == 0 ? nothing : error("numBad(net) should be 0")
    net.numhybrids != 1 ? error("should have 1 hybrid, but net.numhybrids is $(net.numhybrids): $([n.number for n in net.hybrid])") : nothing
end

# Case E bad triangle I
function testCaseE(net::HybridNetwork)
    visited!(net, [istIdentifiable(e) for e in net.edge]);
    n = searchHybridNode(net);
    node = n[1];
    k(node) != 3 ? error("k diff than 3") : nothing
    isVeryBadTriangle(node) ? nothing : error("does not know it is very bad triangle")
    isExtBadTriangle(node) ? error("thinks it is extremely bad triangle") : nothing
    hasVeryBadTriangle(net) ? nothing : error("net does not know it has very bad triangle")
    numBad(net) == 0 ? nothing : error("numBad(net) should be 0")
    net.numhybrids != 1 ? error("should have 1 hybrid, but net.numhybrids is $(net.numhybrids): $([n.number for n in net.hybrid])") : nothing
end


# Case I bad diamond II
function testCaseI(net::HybridNetwork)
    net.numhybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numhybrids)")
    node = net.hybrid[1];
    isBadDiamondII(node) ? nothing : error("does not know it is bad diamond II")
    isBadDiamondI(node) ? error("thinks it is bad diamond I") : nothing
    k(node) == 4 ? nothing : error("k should be 4")
    visited!(net, [istIdentifiable(e) for e in net.edge]);
    edge4 = PhyloNetworks.getIndexEdge(4,net);
    edge1 = PhyloNetworks.getIndexEdge(1,net);
    edge2 = PhyloNetworks.getIndexEdge(2,net);
    edge3 = PhyloNetworks.getIndexEdge(3,net);
    edge9 = PhyloNetworks.getIndexEdge(9,net);
    edge10 = PhyloNetworks.getIndexEdge(10,net);
    edge6 = PhyloNetworks.getIndexEdge(6,net);
    node1 = PhyloNetworks.getIndexNode(-2,net);
    node2 = PhyloNetworks.getIndexNode(-3,net);
    node5 = PhyloNetworks.getIndexNode(-6,net);
    node3 = PhyloNetworks.getIndexNode(3,net);
    (inCycle(net.edge[edge4]) != node.number || inCycle(net.edge[edge9]) != node.number || inCycle(net.edge[edge6]) != node.number || inCycle(net.edge[edge10]) != node.number ) ? error("edges not correctly in cycle") : nothing
    (inCycle(net.node[node1])  != node.number || inCycle(net.node[node2])  != node.number || inCycle(net.node[node5])  != node.number || inCycle(net.node[node3])  != node.number) ? error("nodes 1,5,11,12 not correctly in cycle") : nothing
    (net.edge[edge1].containroot || net.edge[edge2].containroot || net.edge[edge3].containroot) ? error("edges can contain root and shouldn't") : nothing
    (!net.edge[edge4].hybrid || !net.edge[edge4].ismajor) ? error("edge 4 is not hybrid or major") : nothing
    net.edge[edge3].length != 0 ? error("edges should have length 0") : nothing
    istIdentifiable(net.edge[edge3]) ? error("edge9,4 identifiable and should not") : nothing
    (istIdentifiable(net.edge[edge4]) && istIdentifiable(net.edge[edge9]) && istIdentifiable(net.edge[edge6]) && istIdentifiable(net.edge[edge10])) || error("edges that should be identifiable, are not")
    visited(net)[edge4] = false;
    visited(net)[edge9] = false;
    visited(net)[edge6] = false;
    visited(net)[edge10] = false;
    !all([!id for id in visited(net)]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,5,6] ? nothing : error("net.leaf is wrong")
end


# tree example
function testTree(net::HybridNetwork)
    !all([!e.hybrid for e in net.edge]) ? error("some edge is still hybrid") : nothing
    !all([!e.hybrid for e in net.node]) ? error("some node is still hybrid") : nothing
    !all([!hasHybEdge(e) for e in net.node]) ? error("some node has hybrid edge") : nothing
    !all([e.ismajor for e in net.edge]) ? error("some edge is not major") : nothing
    !all([e.containroot for e in net.edge]) ? error("some edge cannot contain root") : nothing
    edge9 = getIndexNumEdge(9,net);
    edge5 = getIndexNumEdge(5,net);
    (!istIdentifiable(net.edge[5]) || !istIdentifiable(net.edge[3])) ? error("edge3,5 not identifiable") : nothing
    visited!(net, [istIdentifiable(e) for e in net.edge]);
    visited(net)[3] = false;
    visited(net)[5] = false;
    !all([!id for id in visited(net)]) ? error("edges not identifiable as identifiable") : nothing
    net.edge[2].length != 1.5 ? error("edge length for 2 is wrong") : nothing
    net.edge[4].length != 0.2 ? error("edge length for 4 is wrong") : nothing
    [n.number for n in net.leaf] == [1,2,4,6,7] ? nothing : error("net.leaf is wrong")
end
