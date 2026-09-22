# SNaQ scores, composite log-likelihood and composite deviance.

"""
    computeSNaQscore!(N, q, ρ=0.0)
    computeSNaQscore!(N, dcf, ρ=0.0)

Computes the composite log-likelihood of network `N` given observed quartet concordance
factor data. The optional `ρ` argument (default 0) is the inheritance correlation parameter
in [0, 1]; `ρ = 0` is independent inheritance, `ρ = 1` is completely dependent.
"""
function computeSNaQscore!(N::HybridNetwork, q::Matrix{Float64}, ρ::Real=0.0)::Float64
    N = deepcopynetwork(N)
    semidirectnetwork!(N)
    qdata, _, params, _, _ = findquartetequations(N)
    loss = computeSNaQscore!(qdata, params, q, ρ)
    SNaQscore!(N, loss)
    return loss
end


function computeSNaQscore!(N::HybridNetwork, dcf::DataCF, ρ::Real=0.0)::Float64
    obsCFs = gatherCFmatrix(dcf)
    eqns, _, parameters, _ = findquartetequations(N);
    loss = computeSNaQscore!(eqns, parameters, obsCFs, ρ)
    SNaQscore!(N, loss)
    for (i, q) in enumerate(dcf.quartet)
        eqn = eqns[i]
        expCF1, expCF2 = computeexpectedCF(eqn, parameters, ρ)
        q.expCF = [expCF1, expCF2, 1.0 - expCF1 - expCF2]
    end
    return loss
end


function computeSNaQscore!(qdata::Vector{QuartetData}, params::Vector{Float64}, q::Matrix{Float64}, ρ::Float64=0.0)::Float64
    α = rhotoalpha(ρ)
    return computelossandgradient!(qdata, params, zeros(length(params)), q, α)
end


"""
    computeSNaQscore!(N, trees, ρ=0.0; propQuartets=1.0, seed=rand(Int))
    computeSNaQscore!(N, lazyq, ρ=0.0; propQuartets=1.0, seed=rand(Int))

Composite log-likelihood of `N` scored directly against gene trees, or against a
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
Computes the composite deviance of the network `net` from the
data in `dcf` with inheritance correlation parameter `ρ`.

Note: if `q.expCF` is populated in each quartet contained in `dcf.quartet`,
then this computation only utilizes `dcf`. Otherwise, expected concordance
factors are computed for `net` first by calling `computeSNaQscore!(net, dcf, ρ)`.
"""
function compositedeviance(net::HybridNetwork, dcf::DataCF, ρ::Float64=0.0)::Float64
    0 <= ρ <= 1 || error("ρ must be in the range [0, 1].")
    if any(q -> q.ngenes <= 0 || ismissing(q.ngenes), dcf.quartet)
        error("At least one quartet in `dcf.quartet` had `q.ngenes` as a value <= 0 or `missing`.")
    end
    if any(q -> length(q.expCF) == 0, dcf.quartet)
        computeSNaQscore!(net, dcf, ρ)
    end

    cdev::Float64 = 0.0
    for q in dcf.quartet
        cdev += q.obsCF[1] == 0.0 ? 0.0 : 2 * q.ngenes * q.obsCF[1] * (log(q.obsCF[1]) - log(q.expCF[1]))
        cdev += q.obsCF[2] == 0.0 ? 0.0 : 2 * q.ngenes * q.obsCF[2] * (log(q.obsCF[2]) - log(q.expCF[2]))
        cdev += q.obsCF[3] == 0.0 ? 0.0 : 2 * q.ngenes * q.obsCF[3] * (log(q.obsCF[3]) - log(q.expCF[3]))
    end
    return cdev
end


"""
Computes the composite log-likelihood of the network `net` from the
data in `dcf` with inheritance correlation parameter `ρ`.

Note: if `q.expCF` is populated in each quartet contained in `dcf.quartet`,
then this computation only utilizes `dcf`. Otherwise, expected concordance
factors are computed for `net` first by calling `computeSNaQscore!(net, dcf, ρ)`.
"""
function compositeloglik(net::HybridNetwork, dcf::DataCF, ρ::Float64=0.0)::Float64
    0 <= ρ <= 1 || error("ρ must be in the range [0, 1].")

    if any(q -> q.ngenes <= 0 || ismissing(q.ngenes), dcf.quartet)
        error("At least one quartet in `dcf.quartet` had `q.ngenes` as a value <= 0 or `missing`.")
    end
    if any(q -> length(q.expCF) == 0, dcf.quartet)
        computeSNaQscore!(net, dcf, ρ)
    end

    cll::Float64 = 0.0
    for q in dcf.quartet
        cll += q.obsCF[1] == 0.0 ? 0.0 : q.ngenes * q.obsCF[1] * log(q.expCF[1])
        cll += q.obsCF[2] == 0.0 ? 0.0 : q.ngenes * q.obsCF[2] * log(q.expCF[2])
        cll += q.obsCF[3] == 0.0 ? 0.0 : q.ngenes * q.obsCF[3] * log(q.expCF[3])
    end
    return cll
end


"""
Computes the composite log-likelihood given the observed concordance factors (`obsCFs`)
and expected concordance factors (`expCFs`) where the number of gene trees used to
calculate the observed CFs in row `i` of `obsCFs` is stored in entry `i` of
`ngt_per_quartet`.
"""
function compositeloglik(obsCFs::Matrix{Float64}, expCFs::Matrix{Float64}, ngt_per_quartet::Vector{Int64})::Float64
    size(obsCFs) == size(expCFs) || error("obsCFs and expCFs have differing sizes.")
    size(obsCFs, 1) == length(ngt_per_quartet) || error("Length of ngt_per_quartet must match the size of the first dimension of obsCFs and expCFs ($(size(obsCFs, 1)) != $(length(ngt_per_quartet)))")
    all(ngt_per_quartet .> 0) || error("All values in ngt_per_quartet must be >0 (min=$(minimum(ngt_per_quartet)))")

    cll::Float64 = 0.0
    for i in eachrow(obsCFs)
        ngt = ngt_per_quartet[i];
        cll += obsCFs[i, 1] == 0.0 ? 0.0 : ngt * obsCFs[i, 1] * log(expCFs[i, 1])
        cll += obsCFs[i, 2] == 0.0 ? 0.0 : ngt * obsCFs[i, 2] * log(expCFs[i, 2])
        cll += obsCFs[i, 3] == 0.0 ? 0.0 : ngt * obsCFs[i, 3] * log(expCFs[i, 3])
    end
    return cll
end


"""
Computes the composite log-likelihood given the observed concordance factors (`obsCFs`)
and expected concordance factors (`expCFs`) assuming all quartets were sampled across
the same number of gene trees `ngt`.
"""
compositeloglik(obsCFs::Matrix{Float64}, expCFs::Matrix{Float64}, ngt::Int64)::Float64 =
    compositeloglik(obsCFs, expCFs, fill(ngt, size(obsCFs, 1)))
