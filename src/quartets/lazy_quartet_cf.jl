# Observed CF matrix that computes its rows on demand.

"""
    LazyQuartetCF(trees, taxa)

Observed quartet concordance factors backed directly by gene trees, standing in for the
materialized matrix [`gatherCFmatrix`](@ref) builds from a `DataCF`.

Row `i` holds the 3 CFs (`12|34, 13|24, 14|23`) of quartet `taxa[unrank4taxa(ntaxa, i)]`,
computed on first access and cached, so no quartet is ever computed unless it is indexed.
That is what lets [`lazysnaq!`](@ref) run on taxon counts where the full `binomial(ntaxa,4)`
matrix would not fit.

`taxa` must be `sort(tiplabels(net))` for the network(s) it is used with -- the same
convention [`findquartetequations!`](@ref) uses -- or ranks will not agree between observed
and expected CFs.

Resolving a row needs only the depths of the 6 pairwise LCAs among its taxa, since the pair
with the deepest LCA is a cherry and the cherry fixes the topology. The constructor
therefore precomputes every pair's LCA depth per gene tree, tree-major so those 6 columns
are read sequentially. That table is `binomial(ntaxa,2) * ntrees` `Int16`s; past
[`MAX_LCADEPTH_BYTES`](@ref) it is skipped and rows are walked per query with
[`observedCF4taxa`](@ref) instead, trading speed for memory linear in the taxon count.
"""
struct LazyQuartetCF <: AbstractMatrix{Float64}
    taxa::Vector{String}                       # canonical, sort(tiplabels(...))-ordered
    ntaxa::Int
    ntrees::Int
    lcadepth::Union{Nothing,Matrix{Int16}}          # (ntrees, npairs); -1 => taxon absent from that tree
    trees::Vector{HybridNetwork}                    # fallback path only (empty otherwise)
    chainmaps::Vector{Dict{String,Vector{Node}}}    # fallback path only (empty otherwise)
    cache::Dict{Int, NTuple{3,Float64}}             # row number => (obsCF12, obsCF13, obsCF14)
    lock::ReentrantLock
end



"""
Largest pairwise-LCA-depth table [`LazyQuartetCF`](@ref) will build before falling back to
walking gene trees per query.
"""
const MAX_LCADEPTH_BYTES::Int = 2 * 1024^3


"""
Cap on the number of quartet rows [`LazyQuartetCF`](@ref) keeps cached; the cache is emptied
once it grows past this.

The cache exists so that materializing `q[rows, :]`, which reads the 3 columns in 3 separate
passes, computes each row once. Runs draw fresh quartets, so there is little to gain from
keeping rows beyond that and much to lose: uncapped, a long multi-run job accumulates
millions of them.
"""
const MAX_CACHED_QUARTETS::Int = 1_000_000


function LazyQuartetCF(trees::Vector{HybridNetwork}, taxa::Vector{String})
    taxa = sort(taxa);
    ntaxa = length(taxa)
    ntaxa >= 4 || error("LazyQuartetCF needs at least 4 taxa (got $(ntaxa)).")
    for tre in trees
        tre.numhybrids == 0 || error("gene tree found that is a network: $(writenewick(tre))")
    end
    ntrees = length(trees)
    npairs = (ntaxa * (ntaxa - 1)) ÷ 2

    if 2 * npairs * ntrees > MAX_LCADEPTH_BYTES
        # Too many taxa to tabulate every pair: fall back to per-query tree walking.
        chainmaps = Vector{Dict{String,Vector{Node}}}(undef, ntrees)
        for (k, tre) in enumerate(trees)
            cmap = Dict{String,Vector{Node}}()
            sizehint!(cmap, length(tre.leaf))
            for l in tre.leaf
                cmap[l.name] = ancestorchain(l)
            end
            chainmaps[k] = cmap
        end
        return LazyQuartetCF(taxa, ntaxa, ntrees, nothing, trees, chainmaps,
                             Dict{Int,NTuple{3,Float64}}(), ReentrantLock())
    end

    taxonindex = Dict{String,Int}(name => i for (i, name) in enumerate(taxa))
    lcadepth = fill(Int16(-1), ntrees, npairs)
    Threads.@threads for t = 1:ntrees
        tre = trees[t]
        buf = Vector{Int}(undef, length(tre.leaf))
        pairlcadepths!(lcadepth, t, tre.node[tre.rooti], 0, buf, 1, taxonindex, ntaxa)
    end

    return LazyQuartetCF(taxa, ntaxa, ntrees, lcadepth, HybridNetwork[],
                         Dict{String,Vector{Node}}[], Dict{Int,NTuple{3,Float64}}(), ReentrantLock())
end


"""
    pairindex(i, j, ntaxa)

Index of the taxon pair `(i,j)`, `i < j`, in [`LazyQuartetCF`](@ref)'s `lcadepth` table.
"""
@inline function pairindex(i::Int, j::Int, ntaxa::Int)::Int
    return ((i - 1) * (2 * ntaxa - i)) ÷ 2 + (j - i)
end


"""
Post-order walk filling row `t` of `D` with the LCA depth of every pair of taxa in one gene
tree. Each pair is written once, at the node where the two lineages first meet. `buf` holds
the taxon indices below the current node, contiguously from `start`; the return value is how
many were written.
"""
function pairlcadepths!(D::Matrix{Int16}, t::Int, node::Node, depth::Int, buf::Vector{Int},
                           start::Int, taxonindex::Dict{String,Int}, ntaxa::Int)::Int
    if node.leaf
        idx = get(taxonindex, node.name, 0)
        idx != 0 || error("gene tree leaf $(node.name) is not in the taxon list.")
        buf[start] = idx
        return 1
    end
    n = 0
    d16 = Int16(depth)
    for ch in getchildren(node)
        m = pairlcadepths!(D, t, ch, depth + 1, buf, start + n, taxonindex, ntaxa)
        # Every pair with one taxon among those already collected for this node and the
        # other among the child just processed has its LCA right here.
        @inbounds for x = start:(start + n - 1)
            bx = buf[x]
            for y = (start + n):(start + n + m - 1)
                by = buf[y]
                a, b = bx < by ? (bx, by) : (by, bx)
                D[t, pairindex(a, b, ntaxa)] = d16
            end
        end
        n += m
    end
    return n
end


Base.size(q::LazyQuartetCF) = (binomial(q.ntaxa, 4), 3)

function Base.show(io::IO, q::LazyQuartetCF)
    print(io, "LazyQuartetCF Object\n")
    print(io, "\tNumber of trees: $(q.ntrees)\n")
    print(io, "\tComputed quartets: $(length(q.cache))/$(binomial(q.ntaxa,4))\n")
end
function Base.show(io::IO, ::MIME"text/plain", q::LazyQuartetCF)
    print(io, "LazyQuartetCF Object\n")
    print(io, "\tNumber of trees: $(q.ntrees)\n")
    print(io, "\tComputed quartets: $(length(q.cache))/$(binomial(q.ntaxa,4))\n")
end

function Base.getindex(q::LazyQuartetCF, i::Int, j::Int)::Float64
    @boundscheck checkbounds(q, i, j)
    cached = lock(q.lock) do
        get(q.cache, i, nothing)
    end
    if cached === nothing
        idx4 = unrank4taxa(q.ntaxa, i)
        D = q.lcadepth
        computed::NTuple{3,Float64} = if D === nothing
            obsCF, _ngenes = observedCF4taxa(q.taxa[idx4], q.trees, q.chainmaps)
            (obsCF[1], obsCF[2], obsCF[3])
        else
            observedCF4taxaidx(idx4[1], idx4[2], idx4[3], idx4[4], D, q.ntaxa, q.ntrees)
        end
        cached = lock(q.lock) do
            length(q.cache) >= MAX_CACHED_QUARTETS && empty!(q.cache)
            get!(q.cache, i, computed)
        end
    end
    return cached[j]
end



"""
    lazyquartetdata(net, trees, propQuartets, seed)
    lazyquartetdata(net, lazyq, propQuartets, seed)

Samples `propQuartets` of `net`'s quartets and materializes just those rows of the observed
CF matrix, returning `(q_idxs, qsub)`. This is what lets the public `lazy` methods score or
fit a network against gene trees without ever building all `binomial(ntaxa,4)` rows.

`net`'s taxa must match the gene trees', and are sorted with `sort` (not
`sort_stringasinteger!`) to match [`findquartetequations!`](@ref)'s convention.
"""
function lazyquartetdata(net::HybridNetwork, lazyq::LazyQuartetCF, propQuartets::Real,
                         seed::Int)::Tuple{Vector{Int},Matrix{Float64}}
    0 < propQuartets <= 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    sort(tiplabels(net)) == lazyq.taxa ||
        error("net's taxa do not match the taxa the observed CFs were built for.")
    q_idxs = sampleqindices(nchoose4taxalength(net), propQuartets, Random.seed!(seed))
    return q_idxs, lazyq[q_idxs, :]
end

function lazyquartetdata(net::HybridNetwork, trees::Vector{HybridNetwork}, propQuartets::Real,
                         seed::Int)::Tuple{Vector{Int},Matrix{Float64}}
    isempty(trees) && error("trees must be non-empty.")
    return lazyquartetdata(net, LazyQuartetCF(trees, sort(tiplabels(net))), propQuartets, seed)
end
