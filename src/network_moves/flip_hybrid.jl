using PhyloNetworks


"""
Randomly selects a hybrid from `N` for which flipping is valid. If no such
hybrids exist in `N`, returns nothing.
"""
function sample_flip_hybrid_parameters(N::HybridNetwork, rng::TaskLocalRNG)::Node
    rperm = randperm(rng, N.numhybrids)
    n0 = deepcopy_network(N)
    for j in rperm
        try
            fliphybrid!(n0, n0.hybrid[j]) === nothing && continue
            return N.hybrid[j]
        catch
            continue
        end
    end
    return nothing
end





