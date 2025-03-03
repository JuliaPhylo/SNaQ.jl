using PhyloNetworks, Test, StatsBase, Random

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

net = reload_labelled_net();
perform_rNNI1!(net, get_nodes(net, "i6", "i4", "i2", "i3")...);
@test writenewick(net) == "(((e,f)i5)#i4,(d,(c,#i4)i6)i3,(a,b)i1)i2;"

net = reload_labelled_net();
perform_rNNI2!(net, get_nodes(net, "i2", "i3", "i6", "i4")...);
@test writenewick(net) == "(((e,f)i5)#i4,(d,(c,#i4)i6)i3,(a,b)i1)i2;"

net = reload_labelled_net();
perform_rNNI3!(net, get_nodes(net, "i6", "f", "i4", "i5")...)
@test writenewick(net) == "((c,#i5)i6,(d,((e)#i5,f)i4)i3,(a,b)i1)i2;"

net = reload_labelled_net();
perform_rNNI4!(net, get_nodes(net, "d", "i6", "i3", "i4")...)
@test writenewick(net) == "((c,#i3)i6,(((e,f)i5,d)i4)#i3,(a,b)i1)i2;"

net = reload_labelled_net();
perform_rNNI4!(net, get_nodes(net, "c", "i3", "i6", "i4")...)
@test !is_valid_rNNI4(get_nodes(net, "i1", "i3", "i2", "i6")...)


##########################################################
# SECOND: VERIFY THAT REVERSIBLE MOVES DO INDEED REVERSE #
##########################################################

# make sure we don't re-use variables
net = nothing
s, t, u, v = nothing, nothing, nothing, nothing
net0 = reload_labelled_net();
net1 = reload_labelled_net();
net2 = reload_labelled_net();
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

net, net0, net1, net2 = nothing, nothing, nothing, nothing
s, t, u, v = nothing, nothing, nothing, nothing
net = reload_labelled_net()

for rNNI_type = 1:4
    for j = 1:100
        prev_newick = writenewick(net)
        valid_stuvs = all_valid_rNNI_nodes(net, rNNI_type)
        if length(valid_stuvs) == 0
            net = reload_labelled_net()
            valid_stuvs = all_valid_rNNI_nodes(net, rNNI_type)
        end
        s, t, u, v = sample(valid_stuvs)

        @test_nowarn perform_rNNI!(net, s, t, u, v, rNNI_type)
        @test_nowarn writenewick(net)
        @test prev_newick != writenewick(net)

        if rNNI_type == 3 && all(getchild(H).leaf for H in net.hybrid)
            # If this is the case, there are 0 valid rNNI(3) moves
            # So, we reset the network
            # (In the case of the current test network, that means
            # we have to reset EVERYTIME...)
            net = reload_labelled_net()
        end
    end
end


net = reload_labelled_net()
for j = 1:1000
    prev_newick = writenewick(net)
    perform_random_rNNI!(net)
    @test_nowarn writenewick(net)
    @test prev_newick != writenewick(net)
end



##############################################################################
# FOURTH: THESE CASES HAVE CAUSED BUGS IN THE PAST - MAKE SURE THEY WORK NOW #
##############################################################################

net = readnewick("(b:1.098,((((f:0.551,(e:0.276)#H1:0.276::0.75)i5:1.146,(c1:2.109,c2:1.557)i6:0.037)i7:0.926,((d1:1.602,d2:1.602)i10:0.086,a:1.687)i12:1.385)i13:1.585,#H1:0.0::0.25)i5554:1.585)i15;")
semidirect_network!(net)
perform_rNNI4!(net, get_nodes(net, "i13", "i5", "i5554", "H1")...)


net = readnewick("(((d1:1.602,d2:1.602)i10:0.086,c2:1.687)i12:1.385,b:4.267,(((c1:2.109,a:1.557)i6:0.037,((e:0.551,f:0.551)i5:0.573)#H1:0.573::0.75)i7:0.463,#H1:0.0::0.25)i-6079:0.463)i13;");
net.isrooted = false;
@test_nowarn all_valid_rNNI1_nodes(net)


net = readnewick("((a:1.146,((c1:2.109,(e:0.551,(f:0.276)#H1:0.5::0.765)i5:0.924)i6:0.005,#H1:0.5::0.235)i8004:0.002)i7:0.0,((d1:1.602,c2:1.602)i10:0.947,d2:1.687)i12:0.0,b:4.267)i13;");
net.isrooted = false;
@test !is_valid_rNNI1(get_nodes(net, "i6", "f", "i8004", "H1")...)


net = readnewick("(a:1.146,(c1:2.109,((e:0.551,(f:0.276)#H1:0.5::0.501)i5:0.425,#H1:0.5::0.499)i3841:0.425)i6:0.007,(c2:0.947,(b:0.947,(d1:1.602,d2:1.602)i10:0.947)i12:0.0)i13:0.0)i7;");
net.isrooted = false;
@test !is_valid_rNNI2(get_nodes(net, "i6", "i5", "i3841", "H1")...)


net = readnewick("((d1:1.602,d2:1.602)i10:0.947,((a:1.146,((c1:2.109,(f:0.551,(e:0.276)#H1:0.5::0.682)i5:1.12)i6:0.0,#H1:0.5::0.318)i-8938:0.007)i7:0.0,c2:1.687)i12:0.0,b:4.267)i13;");
net.isrooted = false;
@test !is_valid_rNNI4(get_nodes(net, "i6", "i5", "i-8938", "H1")...)


net = readnewick("(((b:0.5,(((#H1:0.227::0.491,(c2:0.5,c1:0.875)i3:0.47)i12:0.396,(((d2:0.5,d1:0.75)i7:0.471,((e:0.5,f:1.0)Hi4:0.227)#H1:0.227::0.509)i6216:0.239)#H3:0.239::0.681)i4:0.521,#H3:0.239::0.319)i10076:0.615)i15:0.5,#H2:0.5::0.5)i7554:0.5,(a:0.261)#H2:0.5::0.5)i16;");
net.isrooted = false;
@test_nowarn all_valid_rNNI2_nodes(net)


net = readnewick("((((c1:0.219,c2:1.277)i5355:0.773,((e:0.551,f:0.551)i5:0.48,(#H1:0.5::0.091,(d1:0.65,d2:0.5)i10443:0.052)i8641:0.846)i6:0.126)i7793:1.435,(b:0.441)#H1:0.5::0.909)i7:0.5,a:2.634)i13;");
net.isrooted = false;
s, t, u, v = get_nodes(net, "i13", "i8641", "i7", "H1");
@test !is_valid_rNNI2(s, t, u, v)

