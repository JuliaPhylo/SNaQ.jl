# User-facing entry point taking precomputed concordance factors.

"""
    snaq!(T::HybridNetwork, d::DataCF)
    snaq!(T::HybridNetwork, DataCF(genetrees; lazy=true); propQuartets, propQuartetsFinal)

Estimate the network (or tree) to fit observed quartet concordance factors (CFs)
stored in a DataCF object, using maximum pseudolikelihood.
The search starts from topology `T`,
which can be a tree or a network with no more than `hmax` hybrid nodes.
This function does *not* modify `T`.

If `d` is lazy (see [`DataCF`](@ref)), observed CFs are only computed for the quartets
that are sampled, so `propQuartets` must be less than 1: otherwise every quartet would be
computed, and a `DataCF` with `lazy=false` should be used instead.
Runs then optimize against different samples of quartets, so their resulting networks are
compared on one `propQuartetsFinal` sample shared by every run, without re-optimization,
and `qinfTest` must be false.

Output:

- estimated network in file `.out` (also in `.log`): best network overall and
  list of networks from each individual run.
- the best network and modifications of it, in file `.networks`.
  All networks in this file have the same undirected topology as the best network,
  but have different hybrid/gene flow directions.
  These other networks are reported with their pseudolikelihood scores, because
  non-identifiability issues can cause them to have very similar scores, and
  because SNaQ was shown to estimate the undirected topology accurately but
  not the direction of hybridization in cases of near non-identifiability.
- if any error occurred, file `.err` provides information (seed) to reproduce the error.

There are many optional keyword arguments, including

- `hmax` (default 1): maximum number of hybridizations allowed
- `propQuartets` (default 1): the proportion of observed quartet concordance factors in `d`
  to use when calculating network pseudolikelihoods. Smaller values will lead to faster
  method runtime but may come at the expense of accuracy if lowered too far.
  Must be less than 1 if `d` is lazy.
- `propQuartetsFinal` (default 1): the proportion of quartets in [0, 1], sampled once and
  shared by every run, that the network of each run is re-optimized on at the end
  (see below). If 0, this final re-optimization is skipped.
- `probQR` (default 0): the probability at any given step to use weighted random sampling
  of quartets when deciding where to make topological moves when proposing the next
  candidate network.
- `verbose` (default false): if true, print information about the numerical optimization
- `runs` (default 10): number of independent starting points for the search
- `outgroup` (default none): outgroup taxon to root the estimated topology at the very end
- `filename` (default "snaq"): root name for the output files (`.out`, `.err`). If empty (""),
  files are *not* created, progress log goes to the screen only (standard out).
- `seed` (default 0 to get it from the clock): seed to replicate a given search
- `probST` (default 0.9): probability of perturbing `T` by one NNI move before a run
  starts. With probability 1-probST the run starts from `T` unchanged
  along a tree edge with no hybrid neighbor,
  with a possible modification of one reticulation if `T` has one.
- `updateBL` (default true): If true and if `T` is a tree, the branch lengths in `T`
  are first optimized.

The following optional keyword arguments control when to stop the optimization of branch
lengths and γ's on each individual candidate network. Defaults are in parentheses:

- `ftolRel` (1e-6) and `ftolAbs` (1e-6): relative and absolute differences of
  the network score between the current and proposed parameters,
- `xtolRel` (1e-2) and `xtolAbs` (1e-3): relative and absolute differences
  between the current and proposed parameters.

Greater values will result in a less thorough but faster search.
These parameters are used when evaluating candidate networks only.
The following optional keyword arguments control when to stop proposing new network topologies:

- `Nfail` (50): maximum number of times that new topologies are proposed and rejected (in a row).

Lower values of `Nfail` and greater values of `ftolAbs` would
result in a less thorough but faster search.

At the end, branch lengths and γ's are optimized on the last "best" network of each run
with different and very thorough tolerance parameters:
1e-12 for `ftolRel`, 1e-10 for `ftolAbs`, `xtolRel`, `xtolAbs`.
This uses all quartets if `propQuartetsFinal` is 1 (and only if `propQuartets` < 1),
a `propQuartetsFinal` proportion of quartets if it is in (0, 1), and is skipped if it is 0.

The following optional keyword arguments are used to identify and exclude uninformative quartets.
Uninformative quartets are those with concordance factors sufficiently close to the
expected concordance factors from the star tree (one-third for all topologies). 
Default parameters are in parentheses:
- `qinfTest` (false): if true, then look for uninformative quartets to discard.
- `qtolAbs` (1e-4): tolerance for identifying uninformative concordance factors. Uninformative concordance factors are within `(1/3) ± qtolAbs`.

By default, SNaQ searches for networks under a model of independent inheritance. The following optional
keyword argument controls this model of dependence:

- `ρ` (0.0): inheritance correlation parameter in the range [0, 1]. `ρ = 0` corresponds to
  independent inheritance; `ρ = 1` corresponds to completely dependent inheritance.
  See [Fogg et al. 2023](https://doi.org/10.1093/sysbio/syad030) for further details.

See also: [`fitnumericalparameters!`](@ref) to optimize parameters on a fixed topology,
and [`computeSNaQscore!`](@ref) to get the composite log-likelihood
of a fixed topology with fixed parameters.

References:
  
Claudia Solís-Lemus and Cécile Ané (2016).
Inferring phylogenetic networks with maximum pseudolikelihood under incomplete lineage sorting.
[PLoS Genetics 12(3):e1005896](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)

Kolbow, N, Kong, K, Chafin, T, Justison, J, Ane, C, Solis-Lemus, C (2025).
SNaQ.jl: Improved scalability for phylogenetic network inference.
"""
function snaq!(
  currT0::Union{HybridNetwork, Vector{HybridNetwork}},
  d::DataCF;
  hmax::Int=1,
  Nfail::Int=50,
  ftolRel::Float64=1e-8,
  ftolAbs::Float64=1e-8,
  xtolRel::Float64=1e-8,
  xtolAbs::Float64=1e-8,
  verbose::Bool=false,
  runs::Int=100,
  outgroup::AbstractString="none",
  filename::AbstractString="snaq",
  seed::Int=rand(Int),
  probST::Float64=0.9,
  updateBL::Bool=true,
  probQR::Float64=0.0,
  qtolAbs::Float64=1e-4,
  qinfTest::Bool=false,
  propQuartets::Real=1.0,
  propQuartetsFinal::Real=1.0,
  restrictions::Function=norestrictions,
  ρ::Float64=0.0,
  kwargs...
)
  searchargs = (
      runs=runs,
      maxequivPLs=Nfail,
      verbose=verbose,
      seed=seed,
      probST=probST,
      outgroup=outgroup,
      restrictions=restrictions,
      ftolRel=ftolRel,
      ftolAbs=ftolAbs,
      xtolRel=xtolRel,
      xtolAbs=xtolAbs,
      filename=filename,
      preopt=updateBL,
      qinfTest=qinfTest,
      qtolAbs=qtolAbs,
      probQR=probQR,
      ρ=ρ
  )

  if d.lazy
    propQuartets == 1.0 && error("snaq! with a lazy DataCF requires propQuartets < 1.0, " *
        "otherwise every quartet is computed anyways. Specify propQuartets and " *
        "propQuartetsFinal, or use a DataCF with lazy=false.")
    return lazysnaq!(currT0, d.tree, hmax, propQuartets, propQuartetsFinal; searchargs..., kwargs...)
  end

  bestnet = multisearch(
      currT0,
      d,
      hmax;
      searchargs...,
      propQuartets=propQuartets,
      propQuartetsFinal=propQuartetsFinal,
      kwargs...
  )[1]

  # This call to `fitnumericalparameters!` is only to update the DataCF
  # `d` with the new expected qCFs.
  semidirectnetwork!(bestnet) # for some reason this is being returned with `bestnet.isrooted` as `true`
  fitnumericalparameters!(bestnet, d; maxeval=1)
  for L in bestnet.leaf
    getparentedge(L).length = 0.0
  end
  return bestnet
end
