# SNaQ scores, composite log-likelihood and composite deviance.

"""
    computeSNaQscore!(N, q, ρ=0.0)
    computeSNaQscore!(N, dcf, ρ=0.0)
    computeSNaQscore!(N, dcf, ρ=0.0; propQuartets=0.0, numQuartets=0, seed=42)

Computes the SNaQ score (see [`SNaQscore`](@ref)) of network `N` given observed quartet concordance
factor data. The optional `ρ` argument (default 0) is the inheritance correlation parameter
in [0, 1]; `ρ = 0` is independent inheritance, `ρ = 1` is completely dependent.

If `dcf` is a lazy [`DataCF`](@ref), `N` is scored on a sample of `propQuartets` of all
quartets, or of `numQuartets` quartets, drawn with `seed`, and the CFs of sampled quartets
not yet computed are computed. If neither is given, `N` is scored on the quartets whose CFs
`dcf` has already computed. `propQuartets` and `numQuartets` are not supported otherwise.
"""
function computeSNaQscore!(N::HybridNetwork, q::Matrix{Float64}, ρ::Real=0.0)::Float64
    N = deepcopynetwork(N)
    semidirectnetwork!(N)
    qdata, _, params, _, _ = findquartetequations(N)
    loss = computeSNaQscore!(qdata, params, q, ρ)
    SNaQscore!(N, loss)
    return loss
end


function computeSNaQscore!(N::HybridNetwork, dcf::DataCF, ρ::Real=0.0;
                           propQuartets::Real=0.0, numQuartets::Int=0, seed::Int=42)::Float64
    dcf.lazy && return computeSNaQscorelazy!(N, dcf, Float64(ρ);
                        propQuartets=propQuartets, numQuartets=numQuartets, seed=seed)
    propQuartets == 0 && numQuartets == 0 ||
        error("propQuartets and numQuartets are only supported with a lazy DataCF.")
    obsCFs = gatherCFmatrix(dcf)
    eqns, _, parameters, _ = findquartetequations(N);
    loss = computeSNaQscore!(eqns, parameters, obsCFs, ρ)
    SNaQscore!(N, loss)
    setexpectedCFs!(dcf, N, ρ)
    return loss
end


function computeSNaQscore!(qdata::Vector{QuartetData}, params::Vector{Float64}, q::Matrix{Float64}, ρ::Float64=0.0)::Float64
    α = rhotoalpha(ρ)
    return computelossandgradient!(qdata, params, zeros(length(params)), q, α)
end


"""
    computeSNaQscore!(N, trees, ρ=0.0; propQuartets=1.0, seed=rand(Int))
    computeSNaQscore!(N, lazyq, ρ=0.0; propQuartets=1.0, seed=rand(Int))

SNaQ score (see [`SNaQscore`](@ref)) of `N` scored directly against gene trees, or against a
[`LazyQuartetCF`](@ref) already built from them, without materializing the observed CFs of
every quartet.

`propQuartets` is the proportion of quartets to score on, drawn with `seed`; the default of
1.0 uses them all and so matches `computeSNaQscore!(N, q)`. Pass a smaller value for taxon
counts where all `binomial(ntaxa,4)` quartets will not fit.
"""
function computeSNaQscore!(N::HybridNetwork, trees::Vector{HybridNetwork}, ρ::Real=0.0;
                           propQuartets::Real=1.0, seed::Int=rand(Int))::Float64
    return computeSNaQscore!(N, LazyQuartetCF(trees, sort(tiplabels(N))), ρ;
                             propQuartets=propQuartets, seed=seed)
end

function computeSNaQscore!(N::HybridNetwork, lazyq::LazyQuartetCF, ρ::Real=0.0;
                           propQuartets::Real=1.0, seed::Int=rand(Int))::Float64
    N = deepcopynetwork(N)
    semidirectnetwork!(N)
    q_idxs, qsub = lazyquartetdata(N, lazyq, propQuartets, seed)
    qdata, _, params, _, _ = findquartetequations(N, q_idxs)
    loss = computeSNaQscore!(qdata, params, qsub, Float64(ρ))
    SNaQscore!(N, loss)
    return loss
end


"""
    compositedeviance(net::HybridNetwork, dcf::DataCF, ρ=0.0)
    compositedeviance(obsCFs::Matrix, expCFs::Matrix, ngenes)

Composite deviance of the network `net`, with its current parameters and inheritance
correlation parameter `ρ`, given the observed quartet concordance factors in `dcf`:

    2 * sum over 4-taxon sets q, sum over i=1,2,3 of ngenes_q * oCF_qi * log(oCF_qi / eCF_qi)

where `ngenes_q` is the number of genes for 4-taxon set `q` (field `ngenes` of each quartet
in `dcf`), and terms for which `oCF = 0` contribute 0.
If all 4-taxon sets have the same `ngenes`, then this is `-2 * ngenes * SNaQscore(net)`.

The expected CFs are always calculated from `net`: any CFs stored in `dcf` from an
earlier fit are not used, and `dcf` is not modified.

With matrices, `expCFs` should be in the same order as `obsCFs`, such as obtained from
[`computeexpectedCFmatrix`](@ref) and `SNaQ.gatherCFmatrix`, and `ngenes` can be a single
number, or a vector of the number of genes for each row.
"""
function compositedeviance(net::HybridNetwork, dcf::DataCF, ρ::Real=0.0)::Float64
    obsCFs, expCFs, ngenes = compositeinputs(net, dcf, ρ, "compositedeviance")
    return compositedeviance(obsCFs, expCFs, ngenes)
end

function compositedeviance(obsCFs::AbstractMatrix{<:Real}, expCFs::AbstractMatrix{<:Real},
                           ngenes::Union{Real,AbstractVector{<:Real}})::Float64
    ngenes = checkcompositeinputs(obsCFs, expCFs, ngenes)
    cdev = 0.0
    for i in axes(obsCFs, 1), k in 1:3
        o = obsCFs[i, k]
        o == 0.0 && continue
        cdev += 2 * ngenes[i] * o * (log(o) - log(max(expCFs[i, k], 1e-9)))
    end
    return cdev
end


"""
    compositeloglik(net::HybridNetwork, dcf::DataCF, ρ=0.0)
    compositeloglik(obsCFs::Matrix, expCFs::Matrix, ngenes)

Composite log-likelihood of the network `net`, with its current parameters and inheritance
correlation parameter `ρ`, given the observed quartet concordance factors in `dcf`:

    sum over 4-taxon sets q, sum over i=1,2,3 of ngenes_q * oCF_qi * log(eCF_qi)

where `ngenes_q` is the number of genes for 4-taxon set `q` (field `ngenes` of each quartet
in `dcf`), and terms for which `oCF = 0` contribute 0.

The expected CFs are always calculated from `net`: any CFs stored in `dcf` from an
earlier fit are not used, and `dcf` is not modified.

With matrices, `expCFs` should be in the same order as `obsCFs`, such as obtained from
[`computeexpectedCFmatrix`](@ref) and `SNaQ.gatherCFmatrix`, and `ngenes` can be a single
number, or a vector of the number of genes for each row.
"""
function compositeloglik(net::HybridNetwork, dcf::DataCF, ρ::Real=0.0)::Float64
    obsCFs, expCFs, ngenes = compositeinputs(net, dcf, ρ, "compositeloglik")
    return compositeloglik(obsCFs, expCFs, ngenes)
end

function compositeloglik(obsCFs::AbstractMatrix{<:Real}, expCFs::AbstractMatrix{<:Real},
                         ngenes::Union{Real,AbstractVector{<:Real}})::Float64
    ngenes = checkcompositeinputs(obsCFs, expCFs, ngenes)
    cll = 0.0
    for i in axes(obsCFs, 1), k in 1:3
        o = obsCFs[i, k]
        o == 0.0 && continue
        cll += ngenes[i] * o * log(max(expCFs[i, k], 1e-9))
    end
    return cll
end


"""
Observed CFs, CFs expected from `net`, and number of genes, all in the same order, for
[`compositedeviance`](@ref) and [`compositeloglik`](@ref) (named `fname`) to use.
"""
function compositeinputs(net::HybridNetwork, dcf::DataCF, ρ::Real, fname::AbstractString)
    checknotlazy(dcf, fname)
    0 <= ρ <= 1 || error("ρ must be in the range [0, 1].")
    ngenes = gatherngenes(dcf)
    any(<=(0), ngenes) && error("At least one quartet in `dcf.quartet` has `q.ngenes` <= 0 " *
        "(unknown number of genes): $fname needs the number of genes of every quartet.")
    obsCFs = gatherCFmatrix(dcf)
    expCFs = computeexpectedCFmatrix(net, ρ)
    size(expCFs) == size(obsCFs) || error("the network has $(size(expCFs, 1)) quartets " *
        "but the DataCF has $(size(obsCFs, 1)).")
    return obsCFs, expCFs, ngenes
end

"""
Checks that matrix inputs to [`compositedeviance`](@ref) or [`compositeloglik`](@ref) have
matching sizes, and returns the number of genes as a vector with one value per row.
"""
function checkcompositeinputs(obsCFs::AbstractMatrix, expCFs::AbstractMatrix,
                              ngenes::Union{Real,AbstractVector{<:Real}})::AbstractVector
    size(obsCFs) == size(expCFs) || error("obsCFs and expCFs have differing sizes.")
    size(obsCFs, 2) == 3 || error("obsCFs and expCFs must have 3 columns.")
    ngenes isa Real && (ngenes = fill(ngenes, size(obsCFs, 1)))
    length(ngenes) == size(obsCFs, 1) || error("Length of ngenes must match the number " *
        "of rows of obsCFs and expCFs ($(size(obsCFs, 1)) != $(length(ngenes)))")
    all(>(0), ngenes) || error("All values in ngenes must be > 0 (min = $(minimum(ngenes)))")
    return ngenes
end
