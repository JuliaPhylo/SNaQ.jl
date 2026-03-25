
"""
    snaq!(T::HybridNetwork, d::DataCF)

Estimate the network (or tree) to fit observed quartet concordance factors (CFs)
stored in a DataCF object, using maximum pseudo-likelihood.
The search starts from topology `T`,
which can be a tree or a network with no more than `hmax` hybrid nodes.
This function does *not* modify `T`.

Output:

- estimated network in file `.out` (also in `.log`): best network overall and
  list of networks from each individual run.
- the best network and modifications of it, in file `.networks`.
  All networks in this file have the same undirected topology as the best network,
  but have different hybrid/gene flow directions.
  These other networks are reported with their pseudo-likelihood scores, because
  non-identifiability issues can cause them to have very similar scores, and
  because SNaQ was shown to estimate the undirected topology accurately but
  not the direction of hybridization in cases of near non-identifiability.
- if any error occurred, file `.err` provides information (seed) to reproduce the error.

There are many optional keyword arguments, including

- `hmax` (default 1): maximum number of hybridizations allowed
- `propQuartets` (default 1): the proportion of observed quartet concordance factors in `d`
  to use when calculating network pseudo-likelihoods. Smaller values will lead to faster
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
  are first optimized roughly with [`updateBL!`](@ref) by using the average CF of
  all quartets defining each branch and back-calculating the coalescent units.

The following optional keyword arguments control when to stop the optimization of branch
lengths and γ's on each individual candidate network. Defaults are in parentheses:

- `ftolRel` (1e-6) and `ftolAbs` (1e-6): relative and absolute differences of
  the network score between the current and proposed parameters,
- `xtolRel` (1e-2) and `xtolAbs` (1e-3): relative and absolute differences
  between the current and proposed parameters.

Greater values will result in a less thorough but faster search.
These parameters are used when evaluating candidate networks only.
The following optional keyword arguments control when to stop proposing new network topologies:

- `Nfail` (75): maximum number of times that new topologies are proposed and rejected (in a row).
- `liktolAbs` (1e-6): the proposed network is accepted if its score is better
  than the current score by at least liktolAbs.

Lower values of `Nfail` and greater values of `liktolAbs` and `ftolAbs` would
result in a less thorough but faster search.

At the end, branch lengths and γ's are optimized on the last "best" network
with different and very thorough tolerance parameters:
1e-12 for `ftolRel`, 1e-10 for `ftolAbs`, `xtolRel`, `xtolAbs`.

The following optional keyword arguments are used to identify and exclude uninformative quartets.
Uninformative quartets are those with concordance factors sufficiently close to the
expected concordance factors from the star tree (one-third for all topologies). 
Default parameters are are in parentheses:
- `qinfTest` (false): if true, then look for uninformative quartets to discard.
- `qtolAbs` (1e-4): The tolerance for identifying uninformative concordance factors. Uninformative concordance factors are (1/3)±`qtolAbs`

See also: [`topologymaxQpseudolik!`](@ref) to optimize parameters on a fixed topology,
and [`topologyQpseudolik!`](@ref) to get the deviance (pseudo log-likelihood up to a constant)
of a fixed topology with fixed parameters.

References:
  
Claudia Solís-Lemus and Cécile Ané (2016).
Inferring phylogenetic networks with maximum pseudolikelihood under incomplete lineage sorting.
[PLoS Genetics 12(3):e1005896](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)

Kolbow, N, Kong, K, Chafin, T, Justison, J, Ane, C, Solis-Lemus, C (2025).
SNaQ.jl: Improved scalability for phylogenetic network inference.
"""
function snaq!(
    currT0::HybridNetwork,
    d::DataCF;
    hmax::Integer=1,
    liktolAbs::Float64=1e-8,
    Nfail::Integer=3000,
    ftolRel::Float64=1e-8,
    ftolAbs::Float64=1e-8,
    xtolRel::Float64=1e-8,
    xtolAbs::Float64=1e-8,
    verbose::Bool=false,
    closeN::Bool=true,
    Nmov0::Vector{Int}=Int[],
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
    restrictions::Function=SNaQ.knownidentifiable,
    kwargs...
)

    # kwargs not implemented that should be implemented:
    # - `filename`
    # - `updateBL`
    # - `probQR`
    # - `qinfTest`
    # - `qtolAbs`

    return multisearch(
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
        filename=filename
    )[1]

end