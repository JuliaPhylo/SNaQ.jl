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



