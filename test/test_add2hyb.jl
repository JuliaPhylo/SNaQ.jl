# adding more than one hybrid
# claudia may 2015

tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"
currT0 = readnewicklevel1(tree);
#printedges(currT0)
besttree = deepcopy(currT0);

Random.seed!(16);
successful,hybrid,flag,nocycle,flag2,flag3 = SNaQ.addHybridizationUpdate!(besttree);
@test all([successful, flag, flag2, flag3, !nocycle, hybrid.hybrid])
@test hybrid.number == 11
@test sum(inCycle(e)==11 for e in besttree.edge) == 4
@test k(hybrid) == 4
@test !isBadDiamondI(hybrid)
@test !isBadDiamondII(hybrid)

successful,hybrid,flag,nocycle,flag2,flag3 = SNaQ.addHybridizationUpdate!(besttree); #will add a bad triangle
# seed was chosen such that we tried to add a bad triangle. We should notice.
@test all([!successful, flag, !flag2, flag3, !nocycle])
@test k(hybrid) == 3
@test isVeryBadTriangle(hybrid)
ed = PhyloNetworks.hybridEdges(hybrid)
@test ed[1].ismajor
@test ed[1].gamma > 0.5
@test ed[1].hybrid

# Functions that run stuff w/ the blacklist
deleteHybrid!(besttree.hybrid[2], besttree, true, true)
blacklist(besttree)
@test length(chooseEdgesGamma(besttree, true, besttree.edge)) == 3

tre = readnewicklevel1("(((1,2),3),(4,(5,6)));");
d = readtableCF(joinpath(@__DIR__, "..", "examples", "tableCF.txt"))
Random.seed!(4)
SNaQ.addHybridizationUpdate!(tre)
SNaQ.addHybridizationUpdate!(tre)
deleteHybrid!(tre.hybrid[1], tre, true, true)
@test length(chooseEdgesGamma(tre, true, tre.edge, 0.0, d)) == 3

# did not recognize as bad diamond II
tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
net = readnewicklevel1(tree);
isBadDiamondII(net.node[10]) || error("does not recognize as bad diamond II")
