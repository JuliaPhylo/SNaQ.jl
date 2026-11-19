const EvaluatedTopologies = Dict{Matrix{Int}, Vector{HybridNetwork}}


"""
    rootedattaxon(N)

Copy of `N` rooted at the first available outgroup. Rooted because
`hardwiredclusters` requires rooted networks.
"""
function rootedattaxon(N::HybridNetwork)
    for taxon in sort(tiplabels(N))
        rooted = deepcopynetwork(N)
        try
            rootatnode!(rooted, taxon)
        catch e
            e isa RootMismatch || rethrow()
            continue
        end
        return rooted
    end
    return nothing
end


"""
    topologyin!(N, evaluated)

Checks whether the network `N` is already contained in the dictionary `evaluated`.
Checks by:
1. Comparing hardwired clusters (fastest). Returns false if no hardwired clusters
    in `evaluated` match those of `N`.
2. If a hardwired clusters match is found, compares semi-directed μ-distances
    with `mudistance_semidirected`, which will only return 0 in the case of
    topological isomorphism. Left as a fallback because it is slower.
"""
function topologyin!(N::HybridNetwork, evaluated::EvaluatedTopologies)::Bool
    rooted = rootedattaxon(N)
    isnothing(rooted) && return false
    key = sortslices(hardwiredclusters(rooted, sort(tiplabels(rooted)))[:, 2:end], dims=1)
    sameclusters = get!(evaluated, key, HybridNetwork[])
    any(M -> mudistance_semidirected(M, rooted) == 0, sameclusters) && return true
    push!(sameclusters, rooted)
    return false
end
