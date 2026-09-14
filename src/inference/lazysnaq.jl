# User-facing entry point taking gene trees directly.

"""
    lazysnaq!(currT0, trees, hmax, propQuartets, propQuartetsFinal; kwargs...)

Like [`snaq!`](@ref), but takes gene trees instead of a `DataCF` and never computes the
observed concordance factors of all `binomial(ntaxa,4)` quartets. Only the quartets a search
actually samples are computed, on demand and cached, through a [`LazyQuartetCF`](@ref).

# Arguments
- `currT0`: starting topology, or one per run. If a vector, all must share the same taxa.
- `trees`: gene trees the observed quartet CFs are computed from.
- `hmax`: maximum number of hybridizations allowed.
- `propQuartets`: proportion of quartets each run samples and hill-climbs against.
- `propQuartetsFinal`: proportion of quartets, from one sample shared by every run, used to
  score the runs' resulting networks against each other. Runs optimize against different
  samples, so their own scores are not comparable; this makes the final comparison fair.

# Optional Named Arguments
- `qinfTest` (default `false`): must be `false`. Testing for uninformative quartets needs
  every quartet's observed CF up front; use [`snaq!`](@ref) if you need it.
- `finalopt` (default `false`): must be `false`. Refitting against all quartets would
  materialize every one of them; use [`snaq!`](@ref) if you need it.

All other named arguments match [`snaq!`](@ref).
"""
function lazysnaq!(
  currT0::Union{HybridNetwork, Vector{HybridNetwork}},
  trees::Vector{HybridNetwork},
  hmax::Int,
  propQuartets::Float64,
  propQuartetsFinal::Float64;
  Nfail::Int=100,
  ftolRel::Float64=1e-8,
  ftolAbs::Float64=1e-8,
  xtolRel::Float64=1e-8,
  xtolAbs::Float64=1e-8,
  verbose::Bool=false,
  runs::Int=100,
  outgroup::AbstractString="none",
  filename::AbstractString="lazysnaq",
  seed::Int=rand(Int),
  probST::Float64=0.3,
  updateBL::Bool=true,
  probQR::Float64=0.0,
  qtolAbs::Float64=1e-4,
  qinfTest::Bool=false,
  finalopt::Bool=false,
  restrictions::Function=norestrictions,
  ρ::Float64=0.0,
  kwargs...
)
    qinfTest && error("lazysnaq! does not support qinfTest=true: it needs every " *
        "quartet's observed CF up front. Use snaq! instead.")
    finalopt && error("lazysnaq! does not support finalopt=true: it would materialize " *
        "every quartet. Use snaq! instead.")
    0 < propQuartets <= 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    0 < propQuartetsFinal <= 1 || error("propQuartetsFinal must be in range (0, 1] (propQuartetsFinal = $(propQuartetsFinal))")
    isempty(trees) && error("trees must be non-empty.")

    refnet = currT0 isa HybridNetwork ? currT0 : currT0[1]
    taxa = sort(tiplabels(refnet))  # must match findquartetequations!'s convention
    ntaxa = length(taxa)

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

    ntotal = binomial(ntaxa, 4)
    nusedfinal = propQuartetsFinal == 1.0 ? ntotal : quartetsamplesize(ntotal, propQuartetsFinal)

    lazyq = LazyQuartetCF(trees, taxa)

    _, all_nets = multisearch(
        currT0, lazyq, hmax;
        runs=runs, maxequivPLs=Nfail, verbose=verbose, seed=seed, probST=probST,
        outgroup=outgroup, restrictions=restrictions, ftolRel=ftolRel, ftolAbs=ftolAbs,
        xtolRel=xtolRel, xtolAbs=xtolAbs, propQuartets=propQuartets, filename=filename,
        preopt=updateBL, qinfTest=qinfTest, qtolAbs=qtolAbs, probQR=probQR, ρ=ρ,
        finalopt=finalopt, kwargs...
    )

    finalidxs = sampleqindices(ntotal, propQuartetsFinal, Random.seed!(seed))
    # Materialized once, not per network: every network is scored on the same sample.
    finalq::Matrix{Float64} = lazyq[finalidxs, :]
    finalscores = Vector{Float64}(undef, length(all_nets))
    for (i, net) in enumerate(all_nets)
        N_eqns, _, params, _, _ = findquartetequations(net, finalidxs)
        finalscores[i] = computeSNaQscore!(N_eqns, params, finalq, ρ)
    end

    bestidx = argmax(finalscores)
    bestnet = all_nets[bestidx]
    SNaQscore!(bestnet, finalscores[bestidx])
    logmessage(filename, "lazysnaq!: selected run $bestidx of $runs as the best network " *
        "after fair comparison (SNaQscore = $(round(finalscores[bestidx], digits=5)) on " *
        "$nusedfinal shared quartets).")

    semidirectnetwork!(bestnet) # for some reason this is being returned with `bestnet.isrooted` as `true`
    # Unlike snaq!, no final `fitnumericalparameters!(bestnet, d; maxeval=1)`: that only
    # fills a DataCF's expCF fields, and there is no DataCF here. `bestnet`'s SNaQscore is
    # the one from the shared `propQuartetsFinal` sample.
    for L in bestnet.leaf
        getparentedge(L).length = 0.0
    end

    return bestnet
end
