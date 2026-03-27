using PhyloNetworks


"""
Randomly selects a hybrid from `N` for which flipping is valid. If no such
hybrids exist in `N`, returns nothing.
"""
function samplefliphybridparameters(N::HybridNetwork, rng::TaskLocalRNG)
    rperm = randperm(rng, N.numhybrids)
    for j in rperm
        n0 = deepcopynetwork(N)
        try
            fliphybrid!(n0, n0.hybrid[j]) === nothing && continue
            return Tuple([N.hybrid[j]])
        catch
            continue
        end
    end
    return nothing
end
