using PhyloNetworks, Test

#### Helper functions ####
function reload_labelled_net()
    net = readnewick("((a,b)i1,((c,#H1)i6,(d,((e,f)i5)#H1)i3)i2);")
    semidirect_network!(net)
    net.hybrid[1].name = "i4"   # do it here so that PhyloNetworks doesn't through a warning
    return net
end
get_nodes(net::HybridNetwork, names::String...) = [
    net.node[findfirst(node -> node.name == name, net.node)] for name in names
]

###############################################################
# FIRST: TESTS THAT WHOSE OUTCOMES HAVE BEEN VERIFIED BY HAND #
###############################################################

net = reload_labelled_net()
@test is_valid_rSPR(get_nodes(net, "i1", "i3", "i6", "i2", "i4", "i5")...)
@test !is_valid_rSPR(get_nodes(net, "i1", "i3", "i6", "i2", "i5", "i4")...)
@test !is_valid_rSPR(get_nodes(net, "i3", "i1", "i2", "i6", "i4", "i5")...)

perform_rSPR!(net, get_nodes(net, "i1", "i3", "i6", "i2", "i4", "i5")...)
@test writenewick(net) == "(d,(((e,f)i5,(a,b)i1)i2)#i4,(c,#i4)i6)i3;"

perform_rSPR!(net, get_nodes(net, "i6", "i3", "i2", "i4", "i5", "e")...)
@test writenewick(net) == "(d,((f,(e)#i4)i5,(a,b)i1)i2,(c,#i4)i6)i3;"