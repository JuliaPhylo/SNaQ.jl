using PhyloNetworks, Test
include("../../src/network_moves/identifiability_moves.jl")

bad_net = readnewick("(A,(#H1,(B,(C,(D)#H1))));");
adj_net = readnewick(writenewick(bad_net));
shrinkbaddiamonds!(adj_net)
@test adj_net.numhybrids == 0


bad_net = readnewick("(A,(#H1,((B1,B2),(C,(D)#H1))));");
adj_net = readnewick(writenewick(bad_net));
shrinkbaddiamonds!(adj_net)
@test adj_net.numhybrids == 0


good_net = readnewick("((A1,A2),(#H1,(B,(C,(D)#H1))));");
semidirectnetwork!(good_net);
shrinkbaddiamonds!(good_net)
@test good_net.numhybrids == 1


good_net = readnewick("(A,(#H1,(B,((C1,C2),(D)#H1))));");
semidirectnetwork!(good_net);
shrinkbaddiamonds!(good_net)
@test good_net.numhybrids == 1


good_net = readnewick("(A,(#H1,((B1,B2),(C,((D1,D2))#H1))));");
semidirectnetwork!(good_net);
shrinkbaddiamonds!(good_net)
@test good_net.numhybrids == 1

