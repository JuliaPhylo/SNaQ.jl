# Root-ward paths through a network, and set operations on the edges they cross.

"""
    ancestorchain(node)

`[node, parent(node), grandparent(node), ..., root]`, assuming `node` is in a tree with no
hybrid nodes, so every node has at most 1 parent. Used with [`pathedges`](@ref) to get the
path between two leaves. See [`treelikeancestorchain`](@ref) for the network case.
"""
function ancestorchain(node::Node)::Vector{Node}
    chain = Node[node]
    while true
        pa = getparents(node)
        isempty(pa) && break
        node = pa[1]
        push!(chain, node)
    end
    return chain
end


"""
    treelikeancestorchain(node::Node)::Vector{Node}

Like [`ancestorchain`](@ref), but stops below the first reticulation: the walk ends when
the next node up is a hybrid, or when the current node has anything other than exactly
one parent. Those are precisely the conditions under which
no tree-like path between them exists.
"""
function treelikeancestorchain(node::Node)::Vector{Node}
    chain = Node[node]
    while true
        pa = getparents(node)
        length(pa) == 1 || break
        p = pa[1]
        p.hybrid && break
        node = p
        push!(chain, node)
    end
    return chain
end


"""
    pathedges(chainx::Vector{Node}, chainy::Vector{Node})::Vector{Edge}

The edges on the tree-path connecting the two leaves whose `ancestorchain`s are
`chainx` and `chainy`. Both chains necessarily end at the same node (the tree's root),
so their shared suffix -- found by comparing both chains from the root end inward, an
O(depth) walk with plain `===` comparisons and no hashing or allocation -- is the path
from the root down to the leaves' LCA; everything below that, on each side, is that
leaf's own share of the path.
"""
pathedges(chainx::Vector{Node}, chainy::Vector{Node})::Vector{Edge} =
    pathedges!(Edge[], chainx, chainy)


"""
    pathedges!(dest, chainx, chainy)

[`pathedges`](@ref) into `dest`, which is overwritten and returned, so a caller running
over many quartets can reuse one buffer instead of allocating a path per quartet.
"""
function pathedges!(dest::Vector{Edge}, chainx::Vector{Node}, chainy::Vector{Node})::Vector{Edge}
    ix = length(chainx)
    iy = length(chainy)
    while ix > 1 && iy > 1 && chainx[ix-1] === chainy[iy-1]
        ix -= 1
        iy -= 1
    end
    resize!(dest, (ix - 1) + (iy - 1))
    @inbounds for k in 1:(ix-1)
        dest[k] = getparentedge(chainx[k])
    end
    @inbounds for k in 1:(iy-1)
        dest[(ix-1)+k] = getparentedge(chainy[k])
    end
    return dest
end


"""
    meetingpathedges(chainx, chainy)::Union{Nothing,Vector{Edge}}

Edges of the tree-like path between two leaves whose (possibly reticulation-truncated)
chains are `chainx`/`chainy`, or `nothing` if the chains never meet -- i.e. every path
between the leaves passes through a reticulation. Unlike [`pathedges`](@ref) this cannot
anchor at a shared root, so it searches for the first shared node instead; chains are
only a handful of nodes long, so the scan is cheap.
"""
function meetingpathedges(chainx::Vector{Node}, chainy::Vector{Node})::Union{Nothing,Vector{Edge}}
    edges = Edge[]
    return meetingpathedges!(edges, chainx, chainy) ? edges : nothing
end


"""
    meetingpathedges!(dest, chainx, chainy)::Bool

[`meetingpathedges`](@ref) into `dest`, returning whether the chains meet at all. `dest` is
only overwritten when they do.
"""
function meetingpathedges!(dest::Vector{Edge}, chainx::Vector{Node}, chainy::Vector{Node})::Bool
    @inbounds for ix = 1:length(chainx)
        nx = chainx[ix]
        for iy = 1:length(chainy)
            chainy[iy] === nx || continue
            resize!(dest, (ix - 1) + (iy - 1))
            for k = 1:(ix-1)
                dest[k] = getparentedge(chainx[k])
            end
            for k = 1:(iy-1)
                dest[(ix-1)+k] = getparentedge(chainy[k])
            end
            return true
        end
    end
    return false
end


"""
Edges appearing in both `A` and `B` (in `A`'s order), like `intersect` but without
building the intermediate `Set` that `Base.intersect` allocates.
"""
commonedges(A::Vector{Edge}, B::Vector{Edge})::Vector{Edge} = commonedges!(Edge[], A, B)


"""
[`commonedges`](@ref) into `dest`, which is emptied first and returned.
"""
function commonedges!(dest::Vector{Edge}, A::Vector{Edge}, B::Vector{Edge})::Vector{Edge}
    empty!(dest)
    @inbounds for e in A
        e in B && push!(dest, e)
    end
    return dest
end


"""
Whether the edge sets `A` and `B` (as returned by [`pathedges`](@ref)) share no
edges in common. Equivalent to `isempty(intersect(A, B))`, but never allocates the
intersection itself and short-circuits on the first shared edge.
"""
function aredisjointedges(A::Vector{Edge}, B::Vector{Edge})::Bool
    shortvec, longvec = length(A) <= length(B) ? (A, B) : (B, A)
    for e in shortvec
        e in longvec && return false
    end
    return true
end


"""
Parameter indices of the edges in `A` or `B` (deduplicated), i.e. the `param_map` image
of `union(A, B)`, without allocating the union itself.
"""
function unionedgeparams(A::Vector{Edge}, B::Vector{Edge}, param_map::Dict{Int,Int})::Vector{Int}
    out = Int[]
    @inbounds for e in A
        haskey(param_map, e.number) || continue
        p = param_map[e.number]
        p in out || push!(out, p)
    end
    @inbounds for e in B
        e in A && continue
        haskey(param_map, e.number) || continue
        p = param_map[e.number]
        p in out || push!(out, p)
    end
    return out
end


"""
    QuartetPathScratch()

Edge buffers for one quartet's four connecting paths and their intersection.
[`treelikequartetdata`](@ref) runs once per quartet -- millions of times over a
[`search`](@ref) -- and discards these paths before it returns, so they are reused rather
than reallocated per quartet.
"""
struct QuartetPathScratch
    ab::Vector{Edge}
    cd::Vector{Edge}
    ac::Vector{Edge}
    bd::Vector{Edge}
    common::Vector{Edge}
end


QuartetPathScratch() = QuartetPathScratch(Edge[], Edge[], Edge[], Edge[], Edge[])


"""
    pathscratchpool()

One [`QuartetPathScratch`](@ref) per thread, for the threaded quartet loops. Built by the
caller before it enters the parallel region: threads index this pool, they never grow it.
"""
pathscratchpool()::Vector{QuartetPathScratch} =
    [QuartetPathScratch() for _ = 1:Threads.maxthreadid()]
