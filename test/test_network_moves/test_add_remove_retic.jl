using PhyloNetworks, Test


################################################
# FIRST: WE KNOW THESE MOVES SHOULD BE INVALID #
################################################

net = readnewick("(a,((b,#H1)i2,(c,((d,e)i5)#H1)i4)i3)i1;")
e2 = getparentedge(net.leaf[4])
e1 = getparentedge(getparent(e2))
@test is_valid_add_hybrid(e1, e2, net)
@test !is_valid_add_hybrid(e2, e1, net)


################################################################
# SECOND: add_hybrid! SHOULD BE REVERSIBLE WITH remove_hybrid! #
################################################################

net = readnewick("(a,((b,#H1)i2,(c,((e,d)i5)#H1)i4)i3)i1;");
net0 = deepcopy(net);

e1 = getparentedge(net.leaf[4])
e2 = getparentedge(getparent(e1))

newH = add_hybrid!(e2, e1, net);
remove_hybrid!(newH, net)

# (e,d) gets swapped to (d,e) after the operation, but the Newick is equivalent otherwise
@test replace(writenewick(net), "d" => "e") == replace(writenewick(net0), "d" => "e")