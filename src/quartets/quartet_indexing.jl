# Mapping between quartet ranks and the 4 taxa they stand for, and sampling quartets.

"""
    unrank4taxa(n, rank)

The 4 taxon indices (ascending) of the `rank`-th 4-taxon combination of `n` taxa, in the
same lexicographic order [`incrtaxaidx!`](@ref) produces from `[1,2,3,4]`.

Wraps [`whichQuartet`](@ref), which errors for `n <= 4`, to allow the trivial `n == 4` case.
Both [`findquartetequations!`](@ref) and [`LazyQuartetCF`](@ref) must go through this one
function, or an observed CF could be paired with the wrong quartet's equation.
"""
@inline function unrank4taxa(n::Int, rank::Int)::Vector{Int}
    n == 4 && return [1, 2, 3, 4]
    return whichQuartet(n, rank)
end


"""
Helper function to increment the 4-taxa index within `findquartetequations`.
"""
function incrtaxaidx!(ts::Vector{Int})::Nothing
    ind = findfirst(x -> x>1, diff(ts))
    if ind === nothing ind = 4; end
    ts[ind] += 1
    for j in 1:(ind-1)
        ts[j] = j
    end
end


"""
    nchoose4taxalength(net)

Helper function that calculates how many quartet combinations exist.
"""
@inline function nchoose4taxalength(net::HybridNetwork)::Int
    n = net.numtaxa
    return n * (n-1) * (n-2) * (n-3) ÷ 24
end


"""
    quartetsamplesize(n, p)

How many of `n` quartets a `propQuartets` of `p` asks for. Defined once because the samplers
and [`lazysnaq!`](@ref)'s reported sampling plan must agree.
"""
@inline quartetsamplesize(n::Int, p::Real)::Int = max(ceil(Int, n * p), min(n, 10))


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
    return sort(sample(rng, 1:n, quartetsamplesize(n, p), replace=false))
end


sampleqindices(net::HybridNetwork, p::Real, rng::TaskLocalRNG)::Vector{Int} =
    sampleqindices(nchoose4taxalength(net), p, rng)
