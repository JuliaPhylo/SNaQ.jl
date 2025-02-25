using PhyloNetworks, Test
using PhyloPlots


###############################################################
# FIRST: TESTS THAT WHOSE OUTCOMES HAVE BEEN VERIFIED BY HAND #
###############################################################

function reload_labelled_net()
    net = readnewick("((a,b)i1,((c,#H1)i6,(d,((e,f)i5)#H1)i3)i2);")
    semidirect_network!(net)
    net.hybrid[1].name = "i4"   # do it here so that PhyloNetworks doesn't through a warning
    return net
end
get_nodes(net::HybridNetwork, s::String, t::String, u::String, v::String) = [
    net.node[findfirst(node -> node.name == name, net.node)] for name in [s, t, u, v]
]


net = reload_labelled_net()
s, t, u, v = get_nodes(net, "i6", "i4", "i2", "i3")
perform_rNNI1!(net, s, t, u, v)
@test writenewick(net) == "(((e,f)i5)#i4,(d,(c,#i4)i6)i3,(a,b)i1)i2;"

net = reload_labelled_net()
s, t, u, v = get_nodes(net, "i2", "i3", "i6", "i4")
perform_rNNI2!(net, s, t, u, v)
@test writenewick(net) == "(((e,f)i5)#i4,(d,(c,#i4)i6)i3,(a,b)i1)i2;"

net = reload_labelled_net()
p_edges = [e for e in net.edge if getchild(e).hybrid]
s, t, u, v = get_nodes(net, "i6", "f", "i4", "i5")
perform_rNNI3!(net, s, t, u, v)
@test writenewick(net) == "((c,#i5)i6,(d,((e)#i5,f)i4)i3,(a,b)i1)i2;"

net = reload_labelled_net()
p_edges = [e for e in net.edge if getchild(e).hybrid]
s, t, u, v = get_nodes(net, "d", "i6", "i3", "i4")
perform_rNNI4!(net, s, t, u, v)
@test writenewick(net) == "((c,#i3)i6,(((e,f)i5,d)i4)#i3,(a,b)i1)i2;"


##########################################################
# SECOND: VERIFY THAT REVERSIBLE MOVES DO INDEED REVERSE #
##########################################################

# make sure we don't re-use variables
net = nothing
s, t, u, v = nothing, nothing, nothing, nothing
net0 = reload_labelled_net()
net1 = reload_labelled_net()
net2 = reload_labelled_net()
@test writenewick(net1) == writenewick(net2)

##### rNNI (1) #####
for j = 1:100
    # Perform the move
    perform_rNNI1!(net1, get_nodes(net1, "i6", "i4", "i2", "i3")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI1!(net2, get_nodes(net2, "i6", "i4", "i2", "i3")...)
    @test writenewick(net1) == writenewick(net2)

    # Reverse the move
    perform_rNNI1!(net1, get_nodes(net1, "i4", "i6", "i2", "i3")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI1!(net2, get_nodes(net2, "i4", "i6", "i2", "i3")...)
    @test writenewick(net1) == writenewick(net2)
    @test writenewick(net1) == writenewick(net0)
end

##### rNNI (2) #####
net0, net1, net2 = nothing, nothing, nothing
s, t, u, v = nothing, nothing, nothing, nothing
net0 = reload_labelled_net()
net1 = reload_labelled_net()
net2 = reload_labelled_net()
@test writenewick(net1) == writenewick(net2)

for j = 1:100
    # Perform the move
    perform_rNNI2!(net1, get_nodes(net1, "i2", "i3", "i6", "i4")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI2!(net2, get_nodes(net2, "i2", "i3", "i6", "i4")...)
    @test writenewick(net1) == writenewick(net2)

    # Reverse the move
    perform_rNNI2!(net1, get_nodes(net1, "i3", "i2", "i6", "i4")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI2!(net2, get_nodes(net2, "i3", "i2", "i6", "i4")...)
    @test writenewick(net1) == writenewick(net2)
    @test writenewick(net1) == writenewick(net0)
end

##### rNNI ->(3)->(4) #####
net0, net1, net2 = nothing, nothing, nothing
s, t, u, v = nothing, nothing, nothing, nothing
net0 = reload_labelled_net()
net1 = reload_labelled_net()
net2 = reload_labelled_net()
@test writenewick(net1) == writenewick(net2)

for j = 1:100
    # Perform the move
    perform_rNNI3!(net1, get_nodes(net1, "i6", "f", "i4", "i5")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI3!(net2, get_nodes(net2, "i6", "f", "i4", "i5")...)
    @test writenewick(net1) == writenewick(net2)

    # Reverse the move
    perform_rNNI4!(net1, get_nodes(net1, "f", "i6", "i4", "i5")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI4!(net2, get_nodes(net2, "f", "i6", "i4", "i5")...)
    @test writenewick(net1) == writenewick(net2)
    @test writenewick(net1) == writenewick(net0)
end

##### rNNI ->(4)->(3) #####
net0, net1, net2 = nothing, nothing, nothing
s, t, u, v = nothing, nothing, nothing, nothing
net0 = reload_labelled_net()
net1 = reload_labelled_net()
net2 = reload_labelled_net()
@test writenewick(net1) == writenewick(net2)

for j = 1:100
    # Perform the move
    perform_rNNI4!(net1, get_nodes(net1, "d", "i6", "i3", "i4")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI4!(net2, get_nodes(net2, "d", "i6", "i3", "i4")...)
    @test writenewick(net1) == writenewick(net2)

    # Reverse the move
    perform_rNNI3!(net1, get_nodes(net1, "i6", "d", "i3", "i4")...)
    @test writenewick(net1) != writenewick(net2)
    perform_rNNI3!(net2, get_nodes(net2, "i6", "d", "i3", "i4")...)
    @test writenewick(net1) == writenewick(net2)
    @test writenewick(net1) == writenewick(net0)
end


#######################################################################
# THIRD: RUN A LOT OF NETWORK MOVES AND MAKE SURE THAT NOTHING BREAKS #
#######################################################################

