using PhyloNetworks, Test
import SNaQ:
    isvalidaddhybrid,
    addhybrid!, removehybrid!,
    getparentedge, getparent,
    samplemovereticulateoriginparameters,
    isvalidmovereticulateorigin, movereticulateorigin!



################################################
# FIRST: WE KNOW THESE MOVES SHOULD BE INVALID #
################################################

net = readnewick("(a,((b,#H1)i2,(c,((d,e)i5)#H1)i4)i3)i1;")
e2 = getparentedge(net.leaf[4])
e1 = getparentedge(getparent(e2))
@test isvalidaddhybrid(e1, e2, net)
@test !isvalidaddhybrid(e2, e1, net)


################################################################
# SECOND: addhybrid! SHOULD BE REVERSIBLE WITH removehybrid! #
################################################################

net = readnewick("(a,((b,#H1)i2,(c,((e,d)i5)#H1)i4)i3)i1;");
net0 = deepcopy(net);

e1 = getparentedge(net.leaf[4])
e2 = getparentedge(getparent(e1))

newH = addhybrid!(net, e2, e1);
removehybrid!(net, newH)

# (e,d) gets swapped to (d,e) after the operation, but the Newick is equivalent otherwise
@test replace(writenewick(net), "d" => "e") == replace(writenewick(net0), "d" => "e")






net = readnewick("(((d1:1.602,d2:1.602)i10:0.823,(b:0.463,(a:0.573)#H1:0.5::0.973)i3442:1.191)i13:0.111,((c2:0.757,c1:1.554)i6:0.799,(e:0.776,f:0.776)i5:0.497)i12:0.0,#H1:0.5::0.027)i7;"); net.isrooted = false;
@test_nowarn movereticulateorigin!(net, samplemovereticulateoriginparameters(net, Random.seed!(0))...);

new_origin = getparentedge(net.leaf[3])
H = net.hybrid[1]
@test isvalidmovereticulateorigin(H, new_origin, net)


net = readnewick("(b:1.095,(a:0.268,#H2:0.562::0.5)i3649:0.5,((((e:0.679,f:0.515)i5:0.496,(c2:0.51,c1:0.0)i13:0.8)i6:0.109,(d1:0.65,d2:1.602)H1:0.825)i10737:0.562)#H2:0.562::0.5)i7;");
net.isrooted = false;

samp = samplemovereticulateoriginparameters(net, Random.seed!(0))
@test samp !== nothing


##########################################################################
# removeexcessreticulations! removes hybrids with the smallest minor γ first #
##########################################################################

# H2 has minor γ 0.2, H1 has minor γ 0.3: H2 should be removed first
neth2 = readnewick("(((((D:1,C:1):1)#H2:1::0.8,((A:1,B:1):1)#H1:1::0.7):1,(#H1:1::0.3,E:1):1):1,(#H2:1::0.2,O:1):1);")
SNaQ.semidirectnetwork!(neth2)

# hmax at least the number of hybrids: no change
net = deepcopy(neth2)
SNaQ.removeexcessreticulations!(net, 2)
@test net.numhybrids == 2
@test writenewick(net) == writenewick(neth2)

# hmax = 1: only H2 is removed, H1 and its γs are kept
net = deepcopy(neth2)
SNaQ.removeexcessreticulations!(net, 1)
@test net.numhybrids == 1
H = net.hybrid[1]
taxa = tiplabels(net)
@test sort(taxa[hardwiredcluster(getparentedge(H), taxa)]) == ["A", "B"] # H1 kept, not H2
@test getparentedgeminor(H).gamma ≈ 0.3
@test getparentedge(H).gamma ≈ 0.7

# hmax = 0: both minor edges removed, leaving the major tree
net = deepcopy(neth2)
SNaQ.removeexcessreticulations!(net, 0)
@test net.numhybrids == 0
@test hardwiredclusterdistance(net, majortree(neth2), false) == 0
