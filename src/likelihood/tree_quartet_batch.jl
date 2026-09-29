# Flattened form of the quartets no reticulation divides.

"""
    TreeQuartetBatch

The quartet equations that no reticulation divides, flattened into plain arrays and read
sequentially instead of chasing pointers through `Vector{QuartetData}`.

For such a quartet the three expected CFs are `1 - 2E/3` for the resolved topology and `E/3`
for both others, with `E = exp(-branchsum)`, so `log(E/3) = -branchsum - log 3` and the
observed-CF entropy `Σ q_k log(q_k)` is constant over a fit. Precomputing them leaves one
`exp` and one `log` per quartet per evaluation.

Built once per [`fitnumericalparameters!`](@ref) call by [`treequartetbatch`](@ref).
"""
struct TreeQuartetBatch
    coal::Vector{Int}     # concatenated coal_edges over the tree-like quartets
    off::Vector{Int}      # tree-like quartet i owns coal[off[i]:off[i+1]-1]
    qres::Vector{Float64}   # observed CF of that quartet's RESOLVED topology
    qoth::Vector{Float64}   # observed CFs of the other two topologies, summed
    qlogq::Vector{Float64}  # Σ_k q_k log(q_k), skipping q_k == 0 (constant across evaluations)
    general::Vector{Int}    # indices into the ORIGINAL qdata of quartets a hybrid divides
end


"""
    treequartetbatch(qdata, obsCFs)::TreeQuartetBatch

Flattens the quartets of `qdata` that no reticulation divides into a
[`TreeQuartetBatch`](@ref), recording the rest in its `general` field for the full recursive
evaluation. The split is per quartet, not per network: a reticulation only divides the
quartets whose taxa span it, so most quartets of a large network take the fast path even
when `hmax > 0`.
"""
function treequartetbatch(qdata::Vector{QuartetData}, obsCFs::Matrix{Float64})::TreeQuartetBatch
    nq = length(qdata)
    ntree = 0
    total = 0
    @inbounds for j = 1:nq
        if qdata[j].eqn.division_H == -1
            ntree += 1
            total += length(qdata[j].eqn.coal_edges)
        end
    end

    coal = Vector{Int}(undef, total)
    off = Vector{Int}(undef, ntree + 1)
    qres = Vector{Float64}(undef, ntree)
    qoth = Vector{Float64}(undef, ntree)
    qlogq = Vector{Float64}(undef, ntree)
    general = Vector{Int}(undef, nq - ntree)

    pos = 1
    i = 0
    g = 0
    @inbounds for j = 1:nq
        eqn = qdata[j].eqn
        if eqn.division_H != -1
            g += 1
            general[g] = j
            continue
        end
        i += 1
        off[i] = pos
        for p in eqn.coal_edges
            coal[pos] = Int(p)
            pos += 1
        end
        q1 = obsCFs[j, 1]; q2 = obsCFs[j, 2]; q3 = obsCFs[j, 3]
        wc = eqn.which_coal
        r = wc == 1 ? q1 : (wc == 2 ? q2 : q3)
        qres[i] = r
        qoth[i] = (q1 + q2 + q3) - r
        qlogq[i] = ((q1 > 0) ? q1*log(q1) : 0.0) +
                   ((q2 > 0) ? q2*log(q2) : 0.0) +
                   ((q3 > 0) ? q3*log(q3) : 0.0)
    end
    off[ntree+1] = pos
    return TreeQuartetBatch(coal, off, qres, qoth, qlogq, general)
end


"""
Loss and gradient contribution of a [`TreeQuartetBatch`](@ref)'s hybrid-free quartets,
accumulating into `gradient_storage` (which the caller has already zeroed). Mathematically
identical to the `Vector{QuartetData}` path for these quartets; see
[`TreeQuartetBatch`](@ref) for why it is much cheaper.
"""
@fastmath function treebatchlossandgradient!(batch::TreeQuartetBatch, params::Vector{Float64},
                                             gradient_storage::Vector{Float64})::Float64
    coal = batch.coal; off = batch.off
    qres = batch.qres; qoth = batch.qoth; qlogq = batch.qlogq
    nq = length(qres)
    log3 = log(3.0)
    total = 0.0

    @inbounds for j = 1:nq
        lo = off[j]; hi = off[j+1] - Int(1)

        branchsum = 0.0
        for k = lo:hi
            branchsum += params[coal[k]]
        end
        E = exp(-branchsum)

        res = 1.0 - (2.0/3.0)*E
        oth = E/3.0
        res_c = max(res, 1e-9)
        oth_c = max(oth, 1e-9)

        qr = qres[j]; qo = qoth[j]
        # log(oth_c) is -branchsum - log(3) whenever the clamp did not bite.
        logoth = (oth >= 1e-9) ? (-branchsum - log3) : log(1e-9)
        total += qr*log(res_c) + qo*logoth - qlogq[j]

        gradcoef = qr*((2.0/3.0)*E)/res_c - qo*((1.0/3.0)*E)/oth_c
        for k = lo:hi
            gradient_storage[coal[k]] += gradcoef
        end
    end
    return total
end


"""
Loss and gradient contribution of one quartet no reticulation divides. Its expected CFs
depend on `eqn.coal_edges` only through their sum, so every one of those branches gets the
same gradient contribution.
"""
@fastmath @inline function treequartetlossandgradient!(j::Int, eqn::RecursiveCFEquation, params, q, total_grad, col::Int)::Float64
    coal_edges = eqn.coal_edges
    branchsum::Float64 = 0.0
    @inbounds for p in coal_edges
        branchsum += params[p]
    end
    exp_sum::Float64 = exp(-branchsum)

    wc = eqn.which_coal
    eCF1::Float64 = wc == 1 ? 1 - 2/3*exp_sum : 1/3*exp_sum
    eCF2::Float64 = wc == 2 ? 1 - 2/3*exp_sum : 1/3*exp_sum
    eCF3::Float64 = 1.0 - eCF1 - eCF2

    eCF1 = max(eCF1, 1e-9)
    eCF2 = max(eCF2, 1e-9)
    eCF3 = max(eCF3, 1e-9)

    q1 = q[j, 1]; q2 = q[j, 2]; q3 = q[j, 3]
    total_loss_incr::Float64 =
        ((q1 > 0) ? q1 * log(eCF1 / q1) : 0.0) +
        ((q2 > 0) ? q2 * log(eCF2 / q2) : 0.0) +
        ((q3 > 0) ? q3 * log(eCF3 / q3) : 0.0)

    # d(eCF_k)/d(branch) is `2/3*exp_sum` for the resolved topology and `-1/3*exp_sum`
    # for the other two, identically for every branch in `coal_edges`.
    gradcoef::Float64 =
        q1 * (wc == 1 ? 2/3*exp_sum : -1/3*exp_sum) / eCF1 +
        q2 * (wc == 2 ? 2/3*exp_sum : -1/3*exp_sum) / eCF2 +
        q3 * (wc == 3 ? 2/3*exp_sum : -1/3*exp_sum) / eCF3
    @inbounds for p in coal_edges
        total_grad[p, col] += gradcoef
    end
    return total_loss_incr
end
