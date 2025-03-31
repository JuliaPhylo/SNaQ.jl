"""
Generates a `Dict` mapping all edges and hybrid nodes w/ an associated γ parameter to their
index in the vector that is being optimized.
"""
function generate_optimization_map(net::HybridNetwork)::Tuple{Dict{Union{PN.Edge,PN.Node}, Int}, Dict{Int, Union{PN.Edge,PN.Node}}}
    ngammas = length(net.hybrid)
    nedges = net.numedges - net.numtaxa
    map = Dict{Union{PN.Edge,PN.Node}, Int}()
    j = 1
    for H in net.hybrid
        map[H] = j
        j += 1
    end
    for E in net.edge
        getchild(E).leaf && continue
        # E.hybrid && !E.ismajor && continue
        map[E] = j
        j += 1
    end
    return map, Dict{Int, Union{PN.Edge,PN.Node}}(map[key] => key for key in keys(map))
end
nopt_params(net::HybridNetwork)::Int = length(generate_optimization_map(net)[1])


"""
Function mostly used for debugging, gets the relevant parameters for network `net`.
"""
function get_params(net::HybridNetwork)
    vals = []
    optmap, _ = generate_optimization_map(net)
    for key in keys(optmap)
        if typeof(key) <: PN.Node
            push!(vals, getparentedgeminor(key).gamma)
        else
            push!(vals, key.length)
        end
    end
    return vals
end


"""
Helper function that generates a `Vector{Int}` with `n` indices where each
integer from 1 to `n` has probability `p` of appearing.
"""
function sample_qindices(n::Int, p::Real, rng::TaskLocalRNG)::Vector{Int}
    # We take n*p quartets instead of randomly sampling with
    # probability p so that we don't get any bad edge cases
    np = ceil(Int, n*p)
    np = max(np, min(n, 10))

    bv = falses(n)
    return sort(sample(rng, 1:n, np, replace=false))
end
sample_qindices(net::HybridNetwork, p::Real, rng::TaskLocalRNG)::Vector{Int} =
    sample_qindices(nchoose4taxa_length(net), p, rng)


"""
Helper function that calculates how many quartet combinations exist.
"""
@inline function nchoose4taxa_length(net::HybridNetwork)::Int
    n = net.numtaxa
    return n * (n-1) * (n-2) * (n-3) ÷ 24
end


# TODO: faster delete leaf set?


function deepcopy_network(net::HybridNetwork)::HybridNetwork

    # List of nodes W/O attached edges
    node_map::Dict{Int, Node} = Dict{Int, Node}()
    nodec = Array{Node}(undef, net.numnodes)
    for (j, node) in enumerate(net.node)
        nodec[j] = Node(node.number, node.leaf, node.hybrid)
        nodec[j].name = deepcopy(node.name)
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
    return netc
end


