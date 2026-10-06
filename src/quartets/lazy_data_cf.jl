# Lazy DataCF: observed CFs computed on demand, through a LazyQuartetArray.
# `LazyQuartetArray` and `DataCF` are defined in types.jl.

function Base.getindex(lqa::LazyQuartetArray, i::Int)::Quartet
    lq = lqa.lazyq
    q = lock(() -> get(lqa.quartets, i, nothing), lq.lock)
    q === nothing || return q
    1 <= i <= size(lq, 1) || throw(BoundsError(lqa, i))
    taxa = lq.taxa[unrank4taxa(lq.ntaxa, i)]
    newq = Quartet(i, taxa[1], taxa[2], taxa[3], taxa[4], [lq[i, 1], lq[i, 2], lq[i, 3]])
    return lock(() -> get!(lqa.quartets, i, newq), lq.lock)
end

Base.getindex(lqa::LazyQuartetArray, idxs::AbstractVector{<:Int})::Vector{Quartet} =
    [lqa[i] for i in idxs]

# Only quartets whose CFs are computed count: ranks, not positions, index a LazyQuartetArray.
Base.size(lqa::LazyQuartetArray) = (lock(() -> length(lqa.lazyq.cache), lqa.lazyq.lock),)
Base.eachindex(lqa::LazyQuartetArray) = sort!(lock(() -> collect(keys(lqa.lazyq.cache)), lqa.lazyq.lock))
Base.keys(lqa::LazyQuartetArray) = eachindex(lqa)
@inline Base.haskey(lqa::LazyQuartetArray, i::Int)::Bool = haskey(lqa.lazyq, i)

# The AbstractArray fallback would index ranks 1:length(lqa), computing quartets as it prints.
function Base.show(io::IO, lqa::LazyQuartetArray)
    print(io, "LazyQuartetArray: $(length(lqa)) quartets computed out of $(size(lqa.lazyq, 1))")
end
Base.show(io::IO, ::MIME"text/plain", lqa::LazyQuartetArray) = show(io, lqa)


"""
    LazyDataCF(trees)
    LazyDataCF(filename, trees=HybridNetwork[])

Lazy [`DataCF`](@ref), the same as `DataCF(trees; lazy=true)`.

With `filename`, the observed CFs already computed and saved by `write(filename, d)` for
a lazy DataCF `d` are read back from the file. Without the gene `trees` they came from,
only those CFs are available: any other quartet's CF cannot be computed.
"""
LazyDataCF(trees::Vector{HybridNetwork}) = DataCF(trees; lazy=true)

function LazyDataCF(filename::AbstractString, trees::Vector{HybridNetwork}=HybridNetwork[])
    df = CSV.read(filename, DataFrame)
    rowtaxa(r) = [string(r.taxa1), string(r.taxa2), string(r.taxa3), string(r.taxa4)]
    lazyq = if isempty(trees)
        taxa = sort(unique(reduce(vcat, (rowtaxa(r) for r in eachrow(df)); init=String[])))
        LazyQuartetCF(taxa, length(taxa), 0, nothing, HybridNetwork[],
                      Dict{String,Vector{Node}}[], Dict{Int,NTuple{3,Float64}}(), ReentrantLock())
    else
        LazyQuartetCF(trees)
    end
    for r in eachrow(df)
        # rows are quartet ranks, only valid if the taxa here are those the file was written with
        1 <= r.row <= size(lazyq, 1) && lazyq.taxa[unrank4taxa(lazyq.ntaxa, r.row)] == rowtaxa(r) ||
            error("Quartet $(rowtaxa(r)) in row $(r.row) of $filename does not match the quartet " *
                  "of that rank among taxa $(lazyq.taxa). The file must have been written " *
                  "from a lazy DataCF with the same taxa as " *
                  (isempty(trees) ? "those in the file." : "the gene trees."))
        lazyq.cache[r.row] = (r.CF12_34, r.CF13_24, r.CF14_23)
    end
    return DataCF(LazyQuartetArray(lazyq), trees)
end


"""
    write(filename, d::DataCF)
    write(io, d::DataCF)

Write the observed CFs computed so far by the lazy DataCF `d`, one quartet per row,
to be read back with `LazyDataCF(filename)`.
"""
function Base.write(io::IO, d::DataCF)
    d.lazy || error("write only saves the CFs computed by a lazy DataCF. For a DataCF with " *
        "lazy=false, use CSV.write(filename, tablequartetCF(d)) instead.")
    lq = d.quartet.lazyq
    nb = Base.write(io, "row,taxa1,taxa2,taxa3,taxa4,CF12_34,CF13_24,CF14_23\n")
    for i in eachindex(d.quartet)
        taxa = lq.taxa[unrank4taxa(lq.ntaxa, i)]
        ocfs = lock(() -> lq.cache[i], lq.lock)
        nb += Base.write(io, "$i,$(taxa[1]),$(taxa[2]),$(taxa[3]),$(taxa[4]),$(ocfs[1]),$(ocfs[2]),$(ocfs[3])\n")
    end
    return nb
end

Base.write(filename::AbstractString, d::DataCF) = open(io -> Base.write(io, d), filename, "w")


"""
    computeSNaQscorelazy!(net, d, ρ; propQuartets, numQuartets, seed)

Branch of [`computeSNaQscore!`](@ref) for a lazy DataCF `d`: scores `net` on a sample of
`propQuartets` of all quartets, or of `numQuartets` quartets, drawn with `seed`, computing
the CFs of sampled quartets that are not computed yet. If neither is given (0), scores `net`
on the quartets whose CFs are already computed.
"""
function computeSNaQscorelazy!(net::HybridNetwork, d::DataCF, ρ::Float64;
                               propQuartets::Real, numQuartets::Int, seed::Int)::Float64
    lq = d.quartet.lazyq
    nq = size(lq, 1)
    0 <= propQuartets <= 1 || error("propQuartets must be in the range [0, 1]")
    numQuartets >= 0 || error("numQuartets must not be negative.")
    propQuartets > 0 && numQuartets > 0 && error("Both propQuartets and numQuartets cannot be specified")
    numQuartets <= nq || error("numQuartets ($(numQuartets)) must be less than or equal to " *
        "the number of taxa choose 4 ($(nq))")
    sort(tiplabels(net)) == lq.taxa ||
        error("net's taxa do not match the taxa of the lazy DataCF: $(lq.taxa)")
    propQuartets == 0 && numQuartets == 0 && length(d.quartet) == 0 &&
        error("propQuartets set to 0 and numQuartets set to 0, but the lazy DataCF has no CFs " *
              "computed yet! Either specify propQuartets/numQuartets or load a previous lazy " *
              "DataCF from a file.")

    qidxs::Vector{Int} = if propQuartets > 0
        sample(Xoshiro(seed), 1:nq, Int(floor(nq * propQuartets)), replace=false)
    elseif numQuartets > 0
        sample(Xoshiro(seed), 1:nq, numQuartets, replace=false)
    else
        eachindex(d.quartet)
    end

    # Without gene trees (read from a file without them), only the CFs already computed exist.
    if lq.ntrees == 0 && !all(i -> haskey(lq, i), qidxs)
        error("""
        The lazy DataCF was loaded from a file without gene trees, and you are attempting to compute
        the SNaQ score using quartet CFs that were not already computed. Either: (i) do not specify either
        of propQuartets/numQuartets to use all of the data already computed in the lazy DataCF, or
        (ii) reload the lazy DataCF with the original gene trees.
        """)
    end

    eqns, _, parameters, _ = findquartetequations(net, qidxs)
    loss = computeSNaQscore!(eqns, parameters, lq[qidxs, :], ρ)
    SNaQscore!(net, loss)
    return loss
end
