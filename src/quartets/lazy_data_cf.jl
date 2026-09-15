"""
Basically a helper struct that returns `Quartet` objects on demand
(instead of the 3-tuple oCFs returned by LazyQuartetCF) for the
`LazyDataCF` object.
"""
mutable struct LazyQuartetArray{Quartet} <: AbstractVector{Quartet}
    quartets::Dict{Int,Quartet}
    lazyq::LazyQuartetCF
    LazyQuartetArray(lq::LazyQuartetCF) = new{Quartet}(
		Dict{Int,Quartet}(),
		lq
	)
end

function Base.getindex(lqa::LazyQuartetArray, i::Int)::Quartet
	if haskey(lqa.quartets, i)
        return lqa.quartets[i]
    end

    ocfs = lqa.lazyq[i, 1:3]
    taxa = lqa.lazyq.taxa[unrank4taxa(lqa.lazyq.ntaxa, i)]
	lqa.quartets[i] = Quartet(
        i, taxa[1], taxa[2], taxa[3], taxa[4], ocfs
    )
	return lqa.quartets[i]
end

function Base.getindex(lqa::LazyQuartetArray, idxs::AbstractVector{<:Int})::Vector{Quartet}
	ret = Array{Quartet}(undef, length(idxs))
	for (reti, idx) in enumerate(idxs)
		ret[reti] = lqa[idx]
	end
	return ret
end

function Base.length(lqa::LazyQuartetArray)
	return length(lqa.quartets)
end

function Base.size(lqa::LazyQuartetArray)
	return (length(lqa.quartets),)
end

function Base.eachindex(lqa::LazyQuartetArray)
	return collect(keys(lqa.quartets))
end

struct LazyDataCF
    quartet::LazyQuartetArray

    # Similar parameters to DataCF
	taxa::Vector{String}
	numQuartets::Int
    tree::Vector{HybridNetwork}
    numTrees::Int

    function LazyDataCF(trees::Vector{HybridNetwork})
		taxa = sort(reduce(union, tiplabels(t) for t in trees))
		return new(
			LazyQuartetArray(LazyQuartetCF(trees)),
			taxa,
			binomial(length(taxa), 4),
			trees,
			length(trees)
		)
	end
end

function Base.show(io::IO, ::MIME"text/plain", ldcf::LazyDataCF)
	print(io, 
	"""
	LazyDataCF Object
		$(ldcf.numTrees) gene trees with $(length(ldcf.taxa)) unique taxa
		$(length(ldcf.quartet)) CFs computed out of $(ldcf.numQuartets) total possible
	""")
end

function write(filename::String, ldcf::LazyDataCF)
	open(filename, "w+") do f
		Base.write(f, "row,taxa1,taxa2,taxa3,taxa4,CF12_34,CF13_24,CF14_23\n")
		for i in sort(eachindex(ldcf.quartet))
			taxa = ldcf.taxa[unrank4taxa(length(ldcf.taxa), i)]
			ocfs = ldcf.quartet[i].obsCF
			Base.write(f, "$i,$(taxa[1]),$(taxa[2]),$(taxa[3]),$(taxa[4]),$(ocfs[1]),$(ocfs[2]),$(ocfs[3])\n")
		end
	end
end
