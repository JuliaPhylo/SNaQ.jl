# test to see that deleteHybridizationUpdate undoes all attributes
# prompted by Cecile finding cases when containroot was not updated
# Claudia December 2015

SNaQ.setCHECKNET(true)
SNaQ.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")

@testset "test: delete hybridization" begin
global seed, currT0, besttree, net, successful,hybrid,flag,nocycle,flag2,flag3

seed = 485 # 2738 at v0.14.2

currT0 = readnewicklevel1("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
# warning: the random number generator has a local scope:
# with subsets of tests, the same seed would be re-used over and over.
Random.seed!(seed);
besttree = deepcopy(currT0);

# ===== first hybridization ==========================
successful,_ = addHybridizationUpdate!(besttree);
@test successful
net = deepcopy(besttree);

@test sum(!e.containroot for e in net.edge) == 3 # 3 edges can't contain the root
@test sum(inCycle(e)==11 for e in net.edge) == 6
@test sum(inCycle(n)==11 for n in net.node) == 6
@test length(net.partition) == 6
@test Set([e.number for e in p.edges] for p in net.partition) ==
    Set([[14], [11], [10], [9,7,5,3,1,2,4,6,8], [17], [15]])

# ===== second hybridization ==========================
successful = false
successful,_ = addHybridizationUpdate!(besttree);
@test successful
net = deepcopy(besttree);
@test sum(!e.containroot for e in net.edge) == 8
@test sum(inCycle(e)==11 for e in net.edge) == 6 # edges 12, 13, 16, 18, 19, 20
@test sum(inCycle(n)==11 for n in net.node) == 6
@test sum(inCycle(e)==13 for e in net.edge) == 4 # edges 5, 7, 22, 23
@test sum(inCycle(n)==13 for n in net.node) == 4

# test partition
@test length(net.partition) == 9
@test Set([e.number for e in p.edges] for p in net.partition) ==
    Set([[15], [10], [11], [17], [14], [3,1,2], [21,8,9], [6], [4]])

## # ===== identify containroot for net.node[21]
## net0=deepcopy(net);
## hybrid = net.node[21];
## isBadDiamondI(hybrid)
## nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
## [n.number for n in edgesInCycle] == [23,9,7,5,22] || error("wrong identified cycle")
## [n.number for n in nodesInCycle] == [13,14,-4,-5,-6] || error("wrong identified cycle")
## edgesRoot = identifyContainRoot(net,hybrid);
## [n.number for n in edgesRoot] == [3,1,2] || error("wrong identified contain root")
## [e.containroot for e in edgesRoot] == [false,false,false] || error("edges root should be false contain root")
## edges = hybridEdges(hybrid);
## [e.number for e in edges] == [22,23,3] || error("wrong identified edges")
## undoGammaz!(hybrid,net);
## othermaj = getOtherNode(edges[1],hybrid);
## othermaj.number == -6 || error("wrong othermaj")
## edgesmaj = hybridEdges(othermaj);
## [e.number for e in edgesmaj] == [22,5,4] || error("wrong identified edges")
## edgesmaj[3].containroot || error("wrong edgesmaj[3] contain root")
## undoContainRoot!(edgesRoot);
## [e.containroot for e in edgesRoot] == [true,true,true] || error("edges root should be false contain root")
## deleteHybrid!(hybrid,net,true, false)
## [e.containroot for e in edgesRoot] == [true,true,true] || error("edges root should be false contain root")
## printedges(net)

# ================= delete second hybridization =============================
deleteHybridizationUpdate!(net,net.node[21], false,false);
@test length(net.partition) == 6
@test Set([e.number for e in p.edges] for p in net.partition) ==
    Set([[15], [10], [11], [17], [14], [3,1,2,8,9,6,4,7,5]])
@test sum(!e.containroot for e in net.edge) == 3
@test sum(inCycle(e)==11 for e in net.edge) == 6
@test sum(inCycle(e)==13 for e in net.edge) == 0

# =============== delete first hybridization ===================
deleteHybridizationUpdate!(net,net.node[19]);
checkNet(net)
@test length(net.partition) == 0
#printedges(net)
end
