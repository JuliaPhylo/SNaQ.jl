
"""
    rhotoalpha(ρ)

Convert the inheritance correlation parameter `ρ` ∈ [0, 1] to the internal `α` parameter
used by the likelihood functions. `ρ = 0` (independent inheritance) maps to `α = Inf`;
`ρ = 1` (completely dependent) maps to `α = 0`.

See also: [`alphatorho`](@ref)
"""
rhotoalpha(ρ::Real) = ρ == 0.0 ? Inf : (1.0 - ρ) / ρ

"""
    alphatorho(α)

Convert the internal `α` parameter used by the likelihood functions to the inheritance
correlation parameter `ρ` ∈ [0, 1]. `α = Inf` (independent inheritance) maps to `ρ = 0`;
`α = 0` (completely dependent) maps to `ρ = 1`.

See also: [`rhotoalpha`](@ref)
"""
alphatorho(α::Real) = isinf(α) ? 0.0 : 1.0 / (1.0 + α)


"""
    snaq!(T::HybridNetwork, d::DataCF)

Estimate the network (or tree) to fit observed quartet concordance factors (CFs)
stored in a DataCF object, using maximum pseudolikelihood.
The search starts from topology `T`,
which can be a tree or a network with no more than `hmax` hybrid nodes.
This function does *not* modify `T`.

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
- `probQR` (default 0): the probability at any given step to use weighted random sampling
  of quartets when deciding where to make topological moves when proposing the next
  candidate network.
- `verbose` (default false): if true, print information about the numerical optimization
- `runs` (default 10): number of independent starting points for the search
- `outgroup` (default none): outgroup taxon to root the estimated topology at the very end
- `filename` (default "snaq"): root name for the output files (`.out`, `.err`). If empty (""),
  files are *not* created, progress log goes to the screen only (standard out).
- `seed` (default 0 to get it from the clock): seed to replicate a given search
- `probST` (default 0.3): probability to start from `T` at each given run.
  With problability 1-probST, the search is started from an NNI modification of `T`
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

- `Nfail` (100): maximum number of times that new topologies are proposed and rejected (in a row).

Lower values of `Nfail` and greater values of `ftolAbs` would
result in a less thorough but faster search.

At the end, branch lengths and γ's are optimized on the last "best" network
with different and very thorough tolerance parameters:
1e-12 for `ftolRel`, 1e-10 for `ftolAbs`, `xtolRel`, `xtolAbs`.

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

See also: [`optimize!`](@ref) to optimize parameters on a fixed topology,
and [`computeloss`](@ref) to get the composite log-likelihood
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
  hmax::Integer=1,
  Nfail::Integer=100,
  ftolRel::Float64=1e-8,
  ftolAbs::Float64=1e-8,
  xtolRel::Float64=1e-8,
  xtolAbs::Float64=1e-8,
  verbose::Bool=false,
  runs::Integer=10,
  outgroup::AbstractString="none",
  filename::AbstractString="snaq",
  seed::Integer=rand(Int),
  probST::Float64=0.3,
  updateBL::Bool=true,
  probQR::Float64=0.0,
  qtolAbs::Float64=1e-4,
  qinfTest::Bool=false,
  propQuartets::Float64=1.0,
  restrictions::Function=norestrictions,
  ρ::Float64=0.0,
  kwargs...
)
  bestnet = multisearch(
      currT0,
      d,
      hmax;
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
      propQuartets=propQuartets,
      filename=filename,
      preopt=updateBL,
      qinfTest=qinfTest,
      qtolAbs=qtolAbs,
      probQR=probQR,
      ρ=ρ,
      kwargs...
  )[1]

  # This call to `optimize!` is only to update the DataCF
  # `d` with the new expected qCFs.
  semidirectnetwork!(bestnet) # for some reason this is being returned with `bestnet.isrooted` as `true`
  optimize!(bestnet, d; maxeval=1)
  for L in bestnet.leaf
    getparentedge(L).length = 0.0
  end
  return bestnet
end