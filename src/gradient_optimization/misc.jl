"""
    sampleqindices(n, p, informative, rng)

Helper function that generates a `Vector{Int}` with `n` indices where each
integer from 1 to `n` has probability `p` of appearing. Only samples informative
quartets, so only used when `qinfTest` is `true`.
"""
function sampleqindices(n::Int, p::Real, informative::BitVector, rng::TaskLocalRNG)::Vector{Int}
    valididxs::Vector{Int64} = findall(informative)
    ninform::Int64 = length(valididxs)
    nsamp::Int64 = min(ninform, Int64(ceil(n * p)))
    if length(valididxs) == nsamp
        return valididxs
    end
    return sort(sample(rng, valididxs, nsamp, replace=false))
end
sampleqindices(N::HybridNetwork, p::Real, i::BitVector, rng::TaskLocalRNG) =
    sampleqindices(nchoose4taxalength(N), p, i, rng)


"""
    sampleqindices(n, p, rng)

Helper function that generates a `Vector{Int}` with `n` indices where each
integer from 1 to `n` has probability `p` of appearing.
"""
function sampleqindices(n::Int, p::Real, rng::TaskLocalRNG)::Vector{Int}
    # We take n*p quartets instead of randomly sampling with
    # probability p so that we don't get any bad edge cases
    np = ceil(Int, n*p)
    np = max(np, min(n, 10))

    bv = falses(n)
    return sort(sample(rng, 1:n, np, replace=false))
end
sampleqindices(net::HybridNetwork, p::Real, rng::TaskLocalRNG)::Vector{Int} =
    sampleqindices(nchoose4taxalength(net), p, rng)


"""
    nchoose4taxalength(net)

Helper function that calculates how many quartet combinations exist.
"""
@inline function nchoose4taxalength(net::HybridNetwork)::Int
    n = net.numtaxa
    return n * (n-1) * (n-2) * (n-3) ÷ 24
end


function deepcopynetwork(net::HybridNetwork)::HybridNetwork
    # List of nodes W/O attached edges
    node_map::Dict{Int, Node} = Dict{Int, Node}()
    nodec = Array{Node}(undef, net.numnodes)
    for (j, node) in enumerate(net.node)
        nodec[j] = Node(node.number, node.leaf, node.hybrid)
        nodec[j].name = node.name
        node_map[node.number] = nodec[j]
    end

    # List of edges - also attached the edges to their respective nodes
    edgec = Array{Edge}(undef, net.numedges)
    for (j, e) in enumerate(net.edge)
        n1, n2 = node_map[e.node[1].number], node_map[e.node[2].number]
        edgec[j] = Edge(e.number, e.length, e.hybrid, e.gamma, [n1, n2])
        edgec[j].ischild1 = e.ischild1
        edgec[j].ismajor = e.ismajor
        edgec[j].containroot = e.containroot
        push!(n1.edge, edgec[j])
        push!(n2.edge, edgec[j])
    end

    # List of hybrids
    hybc = Array{Node}(undef, net.numhybrids)
    for (j, hyb) in enumerate(net.hybrid)
        hybc[j] = node_map[hyb.number]
    end

    # List of leaves
    leafc = Array{Node}(undef, length(net.leaf))
    for (j, l) in enumerate(net.leaf)
        leafc[j] = node_map[l.number]
    end

    netc = HybridNetwork()
    netc.numtaxa = net.numtaxa
    netc.numnodes = net.numnodes
    netc.numedges = net.numedges
    netc.node = nodec
    netc.edge = edgec
    netc.leaf = leafc
    netc.rooti = net.rooti
    netc.names = net.names
    netc.hybrid = hybc
    netc.numhybrids = net.numhybrids
    netc.isrooted = net.isrooted
    loglik!(netc, loglik(net))
    return netc
end


