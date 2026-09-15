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

@inline Base.setindex!(lqa::LazyQuartetArray, v::Quartet, i::Int) = (lqa.quartets[i] = v)

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
end

function LazyDataCF(trees::Vector{HybridNetwork})
	taxa = sort(reduce(union, tiplabels(t) for t in trees))
	return LazyDataCF(
		LazyQuartetArray(LazyQuartetCF(trees)),
		taxa,
		binomial(length(taxa), 4),
		trees,
		length(trees)
	)
end

function LazyDataCF(filename::String, genetrees::Vector{HybridNetwork})
	df = CSV.read(filename, DataFrame)

	if length(genetrees) == 0
		taxa = Set{String}()
		for r in eachrow(df)
			push!(taxa, String(r.taxa1))
			push!(taxa, String(r.taxa2))
			push!(taxa, String(r.taxa3))
			push!(taxa, String(r.taxa4))
		end
		taxa = sort(collect(taxa))

		lqcf = LazyQuartetCF(taxa, length(taxa), -1, nothing, HybridNetwork[], [], Dict(), ReentrantLock())
		ldcf = LazyDataCF(
			LazyQuartetArray(lqcf),
			taxa,
			binomial(length(taxa), 4),
			genetrees,
			length(genetrees)
		)
		for r in eachrow(df)
			lqcf.cache[r.row] = (r.CF12_34, r.CF13_24, r.CF14_23)
			ldcf.quartet[r.row] = Quartet(r.row, r.taxa1, r.taxa2, r.taxa3, r.taxa4, [r.CF12_34, r.CF13_24, r.CF14_23])
		end
		return ldcf
	else
		taxa = sort(reduce(union, tiplpabels(t) for t in genetrees))
		ldcf = LazyDataCF(
			LazyQuartetArray(LazyQuartetCF(genetrees)),
			taxa,
			binomial(length(taxa), 4),
			genetrees,
			length(genetrees)
		)
		for r in eachrow(df)
			ldcf.quaret.lazyq.cache[r.row] = [r.CF12_34, r.CF13_24, r.CF14_23]
			ldcf.quartet[r.row] = Quartet(r.row, r.taxa1, r.taxa2, r.taxa3, r.taxa4, [r.CF12_34, r.CF13_24, r.CF14_23])
		end
		return ldcf
	end
end
LazyDataCF(filename::String) = LazyDataCF(filename, HybridNetwork[])

function Base.show(io::IO, ::MIME"text/plain", ldcf::LazyDataCF)
	print(io, 
	"""
	LazyDataCF Object
		$(ldcf.numTrees) gene trees with $(length(ldcf.taxa)) unique taxa
		$(length(ldcf.quartet)) CFs computed out of $(ldcf.numQuartets) total possible
	""")
end

function Base.write(io::IOStream, ldcf::LazyDataCF)
	Base.write(io, "row,taxa1,taxa2,taxa3,taxa4,CF12_34,CF13_24,CF14_23\n")
	for i in sort(eachindex(ldcf.quartet))
		taxa = ldcf.taxa[unrank4taxa(length(ldcf.taxa), i)]
		ocfs = ldcf.quartet[i].obsCF
		Base.write(io, "$i,$(taxa[1]),$(taxa[2]),$(taxa[3]),$(taxa[4]),$(ocfs[1]),$(ocfs[2]),$(ocfs[3])\n")
	end
end

function Base.write(filename::String, ldcf::LazyDataCF)
	open(filename, "w+") do f
		write(f, ldcf)
	end
end

function computeSNaQscore!(net::HybridNetwork, ldcf::LazyDataCF, ρ::Float64=0.0; seed::Int=42, propQuartets::Float64=0.0, numQuartets::Int=0)::Float64
	(propQuartets < 0 || propQuartets > 1) && error("propQuartets must be in the range [0, 1]")
	numQuartets >= 0 || error("numQuartets must not be negative.")
	0 < propQuartets <= 1 && numQuartets > 0 && error("Both propQuartets and numQuartets cannot be specified")
	numQuartets <= binomial(net.numtaxa, 4) || error("numQuartets ($(numQuartets)) must be less than or equal to net.numtaxa choose 4 ($(binomial(net.numtaxa, 4)))")
	propQuartets == 0.0 && numQuartets == 0 && length(ldcf.quartet) == 0 && error("propQuartets set to 0 and numQuartets set to 0, but the LazyDataCF object is empty! Either specify propQuartets/numQuartets or load a previous LazyDataCF object from a file.")

	ntaxa = binomial(net.numtaxa, 4)
	qidxs::Vector{Int} = if propQuartets > 0.0
		sample(Xoshiro(seed), 1:ntaxa, Int(floor(ntaxa * propQuartets)), replace=false)
	elseif numQuartets > 0
		sample(Xoshiro(seed), 1:ntaxa, numQuartets, replace=false)
	else
		collect(eachindex(ldcf.quartet))
	end

	# If the ldcf has 0 trees, it was loaded from a file without the associated gene trees,
	# and we need to check to make sure that all of the selected indices in `qidxs` have
	# already been computed, otherwise inform the user 
	if length(ldcf.tree) == 0
		if !all(i -> haskey(ldcf.quartet.lazyq, i), qidxs)
			error("""
			The LazyDataCF object was loaded from a file without gene trees, and you are attempting to compute
			the SNaQ score using quartet CFs that were not already computed. Either: (i) do not specify either
			of propQuartets/numQuartets to use all of the data already computed in the LazyDataCF object, or
			(ii) reload the LazyDataCF object with the original gene trees.
			""")
		end
	end

	Q = Matrix{Float64}(undef, length(qidxs), 3)
	for i in axes(Q, 1)
		Q[i, :] .= ldcf.quartet[qidxs[i]].obsCF
	end

	eqns, _, parameters, _ = findquartetequations(net, qidxs);
	return computeSNaQscore!(eqns, parameters, Q, ρ)
end
