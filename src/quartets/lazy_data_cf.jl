"""
Basically a helper struct that returns `Quartet` objects on demand
(instead of the 3-tuple oCFs returned by LazyQuartetCF) for the
`LazyDataCF` object.
"""
mutable struct LazyQuartetArray{Quartet} <: AbstractVector{Quartet}
    quartets::Dict{Int,Quartet}
    lazyq::LazyQuartetCF
    LazyQuartetArray(lq::LazyQuartetCF) = new{Quartet}(lq)
end

function Base.getindex(lqa::LazyQuartetArray, i::Int)::Quartet
    if has(quartets, i)
        return quartets[i]
    end

    ocfs = lqa.lazyq[i, 1:3]
    taxa = laq.lazyq.taxa[unrank4taxa(lazyq.ntaxa, i)]
	lqa.quartets[i] = Quartet(
        i, taxa[1], taxa[2], taxa[3], taxa[4], ocfs
    )
	return lqa.quartets[i]
end

struct LazyDataCF
    quartet::LazyQuartetArray

    # All the same parameters as DataCF *except* quartet (b/c we don't
    # always have all quartets) and repSpecies, because
    # this does not support multiple alleles at this point
    numQuartets::Int
    tree::Vector{HybridNetwork}
    numTrees::Int

    LazyDataCF(trees::Vector{HybridNetwork}) = new(LazyQuartetArray(LazyQuartetCF(trees)))
end