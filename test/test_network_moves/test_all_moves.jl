#include("test_includes.jl")

# Here, we repeatedly do the following:
# 1) generate a random starting tree t0 with anywhere from 5-40 tips
# 2) generate a random hmax value that is at most the number of tips in the t0
# 2) perform thousands of network moves on the network
# 3) at each step, ensure that the network topology has changed (by comparing newicks)
# 4) at each step, ensure that every node's name can be found in the networks
#    newick generated from `writenewick` (helps ensure all nodes are traversable)

import SNaQ: deepcopy_network, semidirect_network!,
    generate_move_proposal, apply_move!
using Test, Random, StatsBase

for z = 1:500
    # Generate random data
    rng = Random.seed!(z)
    ntaxa = sample(rng, 5:40)
    nhyb = sample(rng, 0:ntaxa)
    t0 = generate_tree(ntaxa, rand(rng, Int))
    t = deepcopy_network(t0)
    semidirect_network!(t0)
    semidirect_network!(t)
    # Store previous newick to make sure the topology is actually
    # changing everytime HWCD is very expensive so we don't compute it
    prev_newick = writenewick(deepcopy_network(t), round = true)

    # Some relevant info for error messages
    prev_move = :none
    prev_params = ()

    for miter = 1:250
        seed = rand(rng, Int)
        Random.seed!(rng, seed)

        t2 = deepcopy_network(t)
        move = params = -1
        try
            move, params = generate_move_proposal(t2, [], nhyb, rng)
        catch e
            @info """
                \n
                Failed on iteration $(z) on move #$(miter).
                seed = $(seed)
                ntaxa, nhyb = $(ntaxa), $(nhyb)
                Newick prior to move: $(prev_newick)

                Previously applied move: $(prev_move)
                Parameters of previous move:\n\n$(prev_params)
            """
            rethrow(e)
        end

        try
            apply_move!(t2, move, params)

            # Alterations that are performed in the main search function need to be done here as well
            removedegree2nodes!(t2);
            while shrink2cycles!(t2) continue end   # keep shrinking until there is nothing to shrink
            for H in t2.hybrid
                ps = getparents(H)
                length(ps) == 2 || error("Hybrid with $(length(ps)) parents occurred.")
                if ps[1] == ps[2]
                    # This is a 2-cycle that was not caught from the `shrink2cycles!` above, implying
                    # that the 2-cycle goes through the root. This is a bad network, so we skip it
                    continue
                end
            end
            semidirect_network!(t2)
            new_newick = writenewick(deepcopy_network(t2), round = true)
            t = t2

            @test move == :rNNI2 || new_newick != prev_newick
            @test findall(node -> !occursin(node.name, new_newick), t.node) == Int[]
            prev_newick = new_newick

            prev_move = move
            prev_params = params
        catch e
            if move == :rNNI1
                @info """
                    \n
                    Failed on iteration $(z) on move #$(miter).
                    seed = $(seed)
                    ntaxa, nhyb = $(ntaxa), $(nhyb)
                """

                println("""
                \nCode to re-produce the error:

                n = readnewick(\"$(prev_newick)\");
                params = get_nodes(n, \"$(params[1].name)\", \"$(params[2].name)\", \"$(params[3].name)\", \"$(params[4].name)\");
                apply_move!(n, :rNNI1, Tuple(params))

                Newick prior to move: $(prev_newick)
                Proposed move: $(move)
                Proposed parameters:\n\n$(params)
                """)
            else
                @info """
                    \n
                    Failed on iteration $(z) on move #$(miter).
                    seed = $(seed)
                    ntaxa, nhyb = $(ntaxa), $(nhyb)

                    Newick prior to move: $(prev_newick)
                    Proposed move: $(move)
                    Proposed parameters:\n\n$(params)
                """
            end
            rethrow(e)
        end
    end
end


