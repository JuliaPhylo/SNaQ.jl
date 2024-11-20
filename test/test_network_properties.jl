# Not written to fit into the broader `runtests.jl` file, written right now to be run stand-alone

using SNaQ, PhyloNetworks, Graphs, Test
import SNaQ: get_network_level, is_galled_tree, is_galled_network

# Level 3, not galled
net1 = readTopology("(((a,#H1)1,(((b)#H2)#H1,#H3)5)2,((#H2)#H3,c)4)3;)")

# Level 1, galled tree
net2 = readTopology("(((x1,#H1),((x2)#H1,x3)),((x4,#H2),(x3)#H2));")

# Level 2, galled network
net3 = readTopology("((x1,#H1),((((x2)#H1,x3),#H2),((x4)#H2,x5)));")

# Level 4, not galled
net4 = readTopology("((((x1,#H1),(((x2)#H1,x3))#H2),(#H2)#H3),((#H3,#H4),((x4)#H4,x5)));")


@testset "get_network_level" begin
    @test get_network_level(net1) == 3
    @test get_network_level(net2) == 1
    @test get_network_level(net3) == 2
    @test get_network_level(net4) == 4
end

@testset "is_galled_tree" begin
    @test !is_galled_tree(net1)
    @test is_galled_tree(net2)
    @test !is_galled_tree(net3)
    @test !is_galled_tree(net4)
end

@testset "is_galled_network" begin
    @test !is_galled_network(net1)
    @test is_galled_network(net2)
    @test is_galled_network(net3)
    @test !is_galled_network(net4)
end

