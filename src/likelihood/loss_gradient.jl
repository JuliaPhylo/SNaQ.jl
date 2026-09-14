# Composite log-pseudolikelihood and its gradient.

"""
    computelossandgradient!(batch, qdata, params, gradient_storage, q, α)

Loss and gradient using a prebuilt [`TreeQuartetBatch`](@ref): the flattened hybrid-free
quartets go through [`treebatchlossandgradient!`](@ref), and the quartets a reticulation
divides go through the general recursion.
"""
function computelossandgradient!(batch::TreeQuartetBatch, qdata::Vector{QuartetData},
                                 params::Vector{Float64}, gradient_storage::Vector{Float64},
                                 q::Matrix{Float64}, α::Float64)::Float64
    fill!(gradient_storage, 0.0)
    total = treebatchlossandgradient!(batch, params, gradient_storage)

    if !isempty(batch.general)
        np = length(params)
        quartet_grad::Array{Float64,3}, running_grads::Vector{RunningGradient}, total_grad::Matrix{Float64} =
            threadscratch(np)
        fill!(total_grad, 0.0)
        for j in batch.general
            total += quartetlossandgradient!(j, qdata, params, q, α, threadcolumns(),
                quartet_grad, running_grads, total_grad)
        end
        nthreadcols::Int = size(total_grad, 2)
        @inbounds for p = 1:np
            acc = 0.0
            for c = 1:nthreadcols
                acc += total_grad[p, c]
            end
            gradient_storage[p] += acc
        end
    end
    return total
end


"""
Computes expected concordance factors and gradients by recursively passing through `qdata`.
"""
@fastmath function computelossandgradient!(qdata::Vector{QuartetData}, params::Vector{T}, gradient_storage::Vector{T}, q::Matrix{T}, α::T=Inf)::T where T<:Float64

    fill!(gradient_storage, 0.0)
    np::Int = length(params)

    # These annotations must name CONCRETE types: `Array{Float64}` is `Array{Float64,N}
    # where N`, which makes every buffer access in the loops below dynamically dispatched.
    quartet_grad::Array{Float64,3}, running_grads::Vector{RunningGradient}, total_grad::Matrix{Float64} =
        threadscratch(np)
    fill!(total_grad, 0.0)

    nq::Int = length(qdata)
    total::Float64 = 0.0
    if useparallelloop(nq)
        total_loss = Threads.Atomic{Float64}(0.0)
        Threads.@threads for j = 1:nq
            Threads.atomic_add!(total_loss,
                quartetlossandgradient!(j, qdata, params, q, α, threadcolumns(), quartet_grad, running_grads, total_grad))
        end
        total = total_loss[]
    else
        # Serial, so accumulate into a plain local: an atomic add per quartet is pure
        # overhead when there is no other thread.
        for j = 1:nq
            total += quartetlossandgradient!(j, qdata, params, q, α, threadcolumns(), quartet_grad, running_grads, total_grad)
        end
    end

    # Sum the per-thread columns in place; `sum(total_grad, dims=2)[:,1]` would allocate two
    # `np`-sized arrays on every NLopt evaluation.
    nthreadcols::Int = size(total_grad, 2)
    @inbounds for p = 1:np
        acc = 0.0
        for c = 1:nthreadcols
            acc += total_grad[p, c]
        end
        gradient_storage[p] = acc
    end
    return total

end


"""
Loss and gradient contribution of the single quartet `qdata[j]`, accumulated into column
`columns[threadid()]` of `total_grad`.
"""
@fastmath @inline function quartetlossandgradient!(j::Int, qdata, params, q, α, columns, quartet_grad, running_grads, total_grad)::Float64
    eqn = qdata[j].eqn
    if eqn.division_H == -1
        # No reticulation divides this quartet, so its whole equation is one exponential in
        # the branch lengths in `eqn.coal_edges`. With no division `params_seen` stays all
        # false and the running gradient stays all ones, so the general path below reduces
        # to exactly this -- but would sweep every parameter to get there.
        return treequartetlossandgradient!(j, eqn, params, q, total_grad,
                                           columns[Threads.threadid()])
    end

    col = columns[Threads.threadid()]
    # A view, not a copy: copying an `(nparam, 3)` buffer out per quartet per gradient
    # evaluation allocated tens of KB every time.
    iter_grad = @view quartet_grad[:, :, col]
    rg = running_grads[col]
    resetrunninggradient!(rg)

    # Only this quartet's own parameters are written below, so only those need clearing
    # beforehand and reading back afterwards.
    relevant = qdata[j].relevant_params
    @inbounds for p in relevant
        iter_grad[p, 1] = 0.0; iter_grad[p, 2] = 0.0; iter_grad[p, 3] = 0.0
    end

    eCF1::Float64, eCF2::Float64 = computeexpectedCFandgradientrecur!(eqn, params, iter_grad, rg, α)
    eCF3::Float64 = 1.0 - eCF1 - eCF2

    eCF1 = max(eCF1, 1e-9)
    eCF2 = max(eCF2, 1e-9)
    eCF3 = max(eCF3, 1e-9)

    inv_eCF1::Float64 = 1.0 / eCF1
    inv_eCF2::Float64 = 1.0 / eCF2
    inv_eCF3::Float64 = 1.0 / eCF3

    total_loss_incr::Float64 =
        ((q[j, 1] > 0) ? q[j, 1] * log(eCF1 / q[j, 1]) : 0.0) +
        ((q[j, 2] > 0) ? q[j, 2] * log(eCF2 / q[j, 2]) : 0.0) +
        ((q[j, 3] > 0) ? q[j, 3] * log(eCF3 / q[j, 3]) : 0.0)

    q1 = q[j, 1]; q2 = q[j, 2]; q3 = q[j, 3]
    @inbounds for p in relevant
        total_grad[p, col] += q1 * iter_grad[p, 1] * inv_eCF1 +
                              q2 * iter_grad[p, 2] * inv_eCF2 +
                              q3 * iter_grad[p, 3] * inv_eCF3
    end
    return total_loss_incr
end


"""
Debugging function - not to be used internally because it will recompute many things.
"""
function computegradient(net::HybridNetwork, obsCFs::Matrix{Float64}, ρ::Real=0.0)::Vector{Float64}
    α = rhotoalpha(ρ)
    params = gatherparams(net);
    grad = zeros(length(params))
    computelossandgradient!(findquartetequations(net)[1], params, grad, obsCFs, α)
    return grad
end
