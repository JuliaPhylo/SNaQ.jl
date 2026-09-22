# Lazy-DataCF branch of snaq!, taking gene trees directly.

"""
    lazysnaq!(currT0, trees, hmax, propQuartets, propQuartetsFinal; runs, seed, filename, qinfTest, ρ, searchargs...)

Branch of [`snaq!`](@ref) for a lazy [`DataCF`](@ref): never computes the observed
concordance factors of all `binomial(ntaxa,4)` quartets. Only the quartets a search
actually samples are computed, on demand and cached, through a [`LazyQuartetCF`](@ref).

# Arguments
- `currT0`: starting topology, or one per run. If a vector, all must share the same taxa.
- `trees`: gene trees the observed quartet CFs are computed from.
- `hmax`: maximum number of hybridizations allowed.
- `propQuartets`: proportion of quartets each run samples and hill-climbs against.
- `propQuartetsFinal`: proportion of quartets, from one sample shared by every run, used to
  score the runs' resulting networks against each other. Runs optimize against different
  samples, so their own scores are not comparable; this makes the final comparison fair.
  Networks are scored with [`rescoreonsharedquartets!`](@ref), without re-optimization, so
  [`search`](@ref)'s own final re-optimization is skipped.

# Named Arguments
- `qinfTest`: must be `false`. Testing for uninformative quartets needs every quartet's
  observed CF up front.
- `runs`, `seed`, `filename`, `ρ`: as passed to [`multisearch`](@ref), and also used here.
- `searchargs...`: all other arguments passed to [`multisearch`](@ref), as built by `snaq!`.
"""
function lazysnaq!(
  currT0::Union{HybridNetwork, Vector{HybridNetwork}},
  trees::Vector{HybridNetwork},
  hmax::Int,
  propQuartets::Real,
  propQuartetsFinal::Real;
  runs::Int,
  seed::Int,
  filename::AbstractString,
  qinfTest::Bool,
  ρ::Float64,
  searchargs...
)
    qinfTest && error("snaq! with a lazy DataCF does not support qinfTest=true: it needs " *
        "every quartet's observed CF up front. Use a DataCF with lazy=false instead.")
    0 < propQuartets <= 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    0 < propQuartetsFinal <= 1 || error("propQuartetsFinal must be in range (0, 1] (propQuartetsFinal = $(propQuartetsFinal))")
    isempty(trees) && error("trees must be non-empty.")

    if propQuartetsFinal == 1.0
        @warn "propQuartetsFinal is 1.0 (the default value). This may be very time consuming for large datasets."
    end

    refnet = currT0 isa HybridNetwork ? currT0 : currT0[1]
    taxa = sort(tiplabels(refnet))  # must match findquartetequations!'s convention

    if !(currT0 isa HybridNetwork)
        all(n -> sort(tiplabels(n)) == taxa, currT0) ||
            error("All starting topologies in currT0 must share the same set of taxa, " *
                  "so that one LazyQuartetCF's rank <-> taxa mapping is valid for every run.")
    end

    treetaxa = sort(reduce(union, (tiplabels(t) for t in trees); init=String[]))
    extra = setdiff(treetaxa, taxa)
    isempty(extra) || error("Gene trees contain taxa not present in currT0: $(extra).")
    nevershown = setdiff(taxa, treetaxa)
    isempty(nevershown) || error("These taxa in currT0 appear in no gene tree, so every " *
        "quartet involving them would have an observed CF of NaN: $(nevershown).")

    lazyq = LazyQuartetCF(trees, taxa)

    _, all_nets = multisearch(
        currT0, lazyq, hmax;
        runs=runs, seed=seed, filename=filename, qinfTest=qinfTest, ρ=ρ,
        propQuartets=propQuartets, propQuartetsFinal=0.0, searchargs...
    )

    nusedfinal = rescoreonsharedquartets!(all_nets, lazyq, propQuartetsFinal, seed, ρ)
    bestidx = argmax(SNaQscore.(all_nets))
    bestnet = all_nets[bestidx]
    logmessage(filename, "snaq! (lazy DataCF): selected run $bestidx of $runs as the best network " *
        "after fair comparison (SNaQscore = $(round(SNaQscore(bestnet), digits=5)) on " *
        "$nusedfinal shared quartets).")

    semidirectnetwork!(bestnet) # for some reason this is being returned with `bestnet.isrooted` as `true`
    # Unlike the non-lazy branch of snaq!, no final `fitnumericalparameters!(bestnet, d; maxeval=1)`:
    # that only fills a DataCF's expCF fields, and a lazy DataCF has none. `bestnet`'s SNaQscore is
    # the one from the shared `propQuartetsFinal` sample.
    for L in bestnet.leaf
        getparentedge(L).length = 0.0
    end

    return bestnet
end
