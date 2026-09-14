# Equations for quartets no reticulation spans.

"""
    TreeQuartetContext

Per-network scratch that makes [`findquartetequations4taxa`](@ref) cheap to call for the
many quartets of one network: a leaf name => `Node` map and every leaf's
[`ancestorchain`](@ref), so no quartet has to scan `net.leaf` or re-climb to the root.

It works for networks as well as trees. Each leaf's chain stops below the first reticulation
above it, so two leaves' chains meet exactly when a tree-like path joins them; non-treelike
quartets fall back to [`reticulatequartetequations`](@ref) but tree-like quartets get
significant speed and memory savings.
"""
struct TreeQuartetContext
    net::HybridNetwork
    leafmap::Dict{String,Node}
    chains::Dict{String,Vector{Node}}
    hashybrids::Bool
    nodepos::Dict{Int,Int}
    edgepos::Dict{Int,Int}
end


"""
    treequartetcontext(net)::TreeQuartetContext

Builds a [`TreeQuartetContext`](@ref) for `net`.
"""
function treequartetcontext(net::HybridNetwork)::TreeQuartetContext
    leafmap = Dict{String,Node}()
    chains = Dict{String,Vector{Node}}()
    sizehint!(leafmap, length(net.leaf))
    sizehint!(chains, length(net.leaf))
    hashybrids = net.numhybrids != 0
    for l in net.leaf
        leafmap[l.name] = l
        chains[l.name] = hashybrids ? treelikeancestorchain(l) : ancestorchain(l)
    end
    nodepos = Dict{Int,Int}()
    edgepos = Dict{Int,Int}()
    if hashybrids     # only the reticulate route needs these
        sizehint!(nodepos, length(net.node))
        sizehint!(edgepos, length(net.edge))
        for (j, n) in enumerate(net.node); nodepos[n.number] = Int(j); end
        for (j, e) in enumerate(net.edge); edgepos[e.number] = Int(j); end
    end
    return TreeQuartetContext(net, leafmap, chains, hashybrids, nodepos, edgepos)
end


"""
[`trytreelikequartet`](@ref) for a hybrid-free network, using a precomputed
[`TreeQuartetContext`](@ref). With no hybrids anywhere, every path between two leaves is
tree-like, so this always succeeds (it never returns `nothing`), and the paths come from
the context's ancestor chains instead of being re-climbed per quartet.
"""
function trytreelikequartet(ctx::TreeQuartetContext, taxa::AbstractVector{String}, param_map::Dict{Int,Int},
                            scratch::QuartetPathScratch=QuartetPathScratch())::Union{QuartetData,Nothing}
    chains = ctx.chains
    return treelikequartetdata(chains[taxa[1]], chains[taxa[2]], chains[taxa[3]], chains[taxa[4]],
                                taxa, param_map, ctx.hashybrids, scratch)
end


"""
The shared body of both [`trytreelikequartet`](@ref) methods, given the four leaves'
ancestor chains. Returns `nothing` when some pair is joined only through a reticulation,
i.e. when the quartet is not tree-like.
"""
function treelikequartetdata(ca::Vector{Node}, cb::Vector{Node}, cc::Vector{Node}, cd::Vector{Node},
                              taxa::AbstractVector{String}, param_map::Dict{Int,Int},
                              hashybrids::Bool,
                              scratch::QuartetPathScratch=QuartetPathScratch())::Union{QuartetData,Nothing}
    # NOTE: the quartet's internal branch is NOT "each cherry's LCA up to the four-taxon
    # MRCA": that picks up extra edges when the root lies inside the quartet's span.
    # Intersecting the two connecting paths, as below, is correct in every rooting.
    path_ab::Vector{Edge} = scratch.ab
    path_cd::Vector{Edge} = scratch.cd
    path_ac::Vector{Edge} = scratch.ac
    path_bd::Vector{Edge} = scratch.bd
    if hashybrids
        # Any pair joined only through a reticulation means this quartet is not tree-like.
        meetingpathedges!(path_ab, ca, cb) || return nothing
        meetingpathedges!(path_cd, cc, cd) || return nothing
        meetingpathedges!(path_ac, ca, cc) || return nothing
        meetingpathedges!(path_bd, cb, cd) || return nothing
    else
        pathedges!(path_ab, ca, cb)
        pathedges!(path_cd, cc, cd)
        pathedges!(path_ac, ca, cc)
        pathedges!(path_bd, cb, cd)
    end

    names = [taxa[1], taxa[2], taxa[3], taxa[4]]
    if aredisjointedges(path_ab, path_cd)
        i_acbd = commonedges!(scratch.common, path_ac, path_bd)
        return QuartetData(
            RecursiveCFEquation(true, [param_map[e.number] for e in i_acbd], 1, -1, EMPTY_EQN_VEC, length(param_map)),
            unionedgeparams(path_ac, path_bd, param_map),
            names
        )
    elseif aredisjointedges(path_ac, path_bd)
        i_abcd = commonedges!(scratch.common, path_ab, path_cd)
        return QuartetData(
            RecursiveCFEquation(true, [param_map[e.number] for e in i_abcd], 2, -1, EMPTY_EQN_VEC, length(param_map)),
            unionedgeparams(path_ab, path_cd, param_map),
            names
        )
    else
        i_abcd = commonedges!(scratch.common, path_ab, path_cd)
        return QuartetData(
            RecursiveCFEquation(true, [param_map[e.number] for e in i_abcd], 3, -1, EMPTY_EQN_VEC, length(param_map)),
            unionedgeparams(path_ab, path_cd, param_map),
            names
        )
    end
end


"""
Uses simple path-finding operations to try and find a tree-like quartet
relationship between the 4 taxa in `taxa`. On successful finding of this
quartet, the corresponding `QuartetData` object is returned. If a hybrid
is encountered along a given path in this operation, `nothing` is
returned instead.
"""
function trytreelikequartet(net::HybridNetwork, taxa::AbstractVector{String}, param_map::Dict{Int,Int})::Union{QuartetData,Nothing}
    # Climb each of the 4 leaves once, rather than each of the 4 pairs separately: the
    # reticulation recursion calls this at every leaf of an equation tree.
    hashybrids = net.numhybrids != 0
    chain(name) = begin
        i = findfirst(l -> l.name == name, net.leaf)
        i === nothing && error("trytreelikequartet: taxon $name is not a leaf of the network.")
        l = net.leaf[i]
        hashybrids ? treelikeancestorchain(l) : ancestorchain(l)
    end
    return treelikequartetdata(chain(taxa[1]), chain(taxa[2]), chain(taxa[3]), chain(taxa[4]),
                                taxa, param_map, hashybrids)
end
