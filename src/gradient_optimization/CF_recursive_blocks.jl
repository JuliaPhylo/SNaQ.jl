


"""
The key struct used in computing -log pseudo-likelihood and gradients
when optimizing branch lengths.
"""
mutable struct RecursiveCFEquation
    can_coalesce_here::Bool
    coal_edges::Vector{Int}
    which_coal::Int     # 0 = NA, 1 = ab|cd, 2 = ac|bd, 3 = ad|bc
    division_H::Int     # hybrid index in `net.hybrid`; -1 means not set
    divisions::Vector{RecursiveCFEquation}
    coal_mask::BitVector

    function RecursiveCFEquation(can_coal::Bool, coal_Es::Vector{Int}, which_c::Int, dH::Int, d::Vector{RecursiveCFEquation}, nparam::Int)
        mask = falses(nparam)
        @inbounds @simd for pidx in coal_Es
            mask[pidx] = true
        end
        new(can_coal, coal_Es, which_c, dH, d, mask)
    end
end


"""
A struct that contains:
1. The initial `RecrusiveCFEquation` struct from which the loss & gradient can be calculated
2. A list of "internal" parameters (stored as indexed from 1 to `k` where `k` is the total number of
    parameters optimized during branch length optimization in `optimize_bls!`) that are relevant to
    this quarnet. This INCLUDES edges that DO NOT contribute to the quarnet's expected CF, but that
    ARE internal edges w/in the network as a whole and DO **inscribe** the quarnet in the network.
    I.e., if one of these edges was removed or re-directed, this quarnet's eCF would change.
3. The set of taxa that these equations relate to.
"""
mutable struct QuartetData
    eqn::RecursiveCFEquation
    relevant_params::Vector{Int}
    q_taxa::SizedVector{4,String}
end
contains_parameter(qdata::QuartetData, param_idxs::Vector{Int})::Bool = any(
    obj_idx -> obj_idx in qdata.relevant_params, param_idxs
)


function compute_eCF(qdata::QuartetData, params::Vector{Float64}, α::Float64)
    return compute_eCF_and_gradient_recur!(qdata.eqn, params, zeros(length(params), 3), falses(length(params)), α, ones(length(params), 3))
end


"""
Helper function (primarily for the `QuartetNetworkGoodnessFit.jl` package) used
to compute the expected CFs of the quartet consisting of `taxa` in `net`.
"""
function compute_eCF_4taxa(net::HybridNetwork, taxa::AbstractVector{<:AbstractString}, α::Real)::Tuple{Float64,Float64,Float64}
    # Definitely slightly inefficient to do this for each quartet, but shouldn't be a big deal.
    param_map, params = gather_optimization_info(net)[[2,4]]

    qdata = SNaQ.find_quartet_equations_4taxa(net, taxa, param_map)
    eCF1, eCF2 = SNaQ.compute_eCF(qdata, params, α)
    
    return eCF1, eCF2, 1-eCF1-eCF2
end


"""
Computes the loss (-log pseudo-likelihood) of network `N` given observed quartet concordance
factor data `q` under Dirichlet parameter `α`.
"""
function compute_loss(N::HybridNetwork, q::Matrix{Float64}, α::Real=Inf)::Float64
    N = deepcopy_network(N)
    semidirect_network!(N)
    qdata, _, params, _, _ = find_quartet_equations(N)
    return compute_loss(qdata, params, q, α)
end
function compute_loss(qdata::Vector{QuartetData}, params::Vector{Float64}, q::Matrix{Float64}, α::Float64=Inf)::Float64
    return compute_loss_and_gradient!(qdata, params, zeros(length(params)), q, α)
end


# Global thread buffers
const THREAD_ITER_GRAD_BUFFER::Dict{Int, Array{Float64, 3}} = Dict{Int, Array{Float64, 3}}()
const THREAD_BV_BUFFER::Dict{Int, BitMatrix} = Dict{Int, BitMatrix}()
const THREAD_RUNNING_GRAD_BUFFER::Dict{Int, Array{Float64, 3}} = Dict{Int, Array{Float64, 3}}()
const THREAD_LOCAL_GRAD_BUFFER::Dict{Int, Array{Float64}} = Dict{Int, Array{Float64}}()

function get_or_create_buffers(params_len::Int)
    if !haskey(THREAD_LOCAL_GRAD_BUFFER, params_len)
        THREAD_ITER_GRAD_BUFFER[params_len] = Array{Float64}(undef, params_len, 3, Threads.nthreads()+1)
        THREAD_BV_BUFFER[params_len] = BitArray(undef, params_len, Threads.nthreads()+1)
        THREAD_RUNNING_GRAD_BUFFER[params_len] = Array{Float64}(undef, params_len, 3, Threads.nthreads()+1)
        THREAD_LOCAL_GRAD_BUFFER[params_len] = Array{Float64}(undef, params_len, Threads.nthreads()+1)
    end
    return THREAD_ITER_GRAD_BUFFER[params_len], THREAD_BV_BUFFER[params_len], THREAD_RUNNING_GRAD_BUFFER[params_len], THREAD_LOCAL_GRAD_BUFFER[params_len]
end


"""
Computes expected concordance factors and gradients by recursively passing through `qdata`.
"""
@fastmath function compute_loss_and_gradient!(qdata::Vector{QuartetData}, params::Vector{T}, gradient_storage::Vector{T}, q::Matrix{T}, α::T=Inf)::T where T<:Float64

    thread_lock::ReentrantLock = ReentrantLock()
    fill!(gradient_storage, 0.0)
    total_loss = Threads.Atomic{Float64}(0.0)
    np::Int = length(params)

    iter_grad_buffer::Array{Float64}, bv_buffer::BitMatrix, running_grad_buffer::Array{Float64}, local_grad_buffer::Array{Float64} =
        get_or_create_buffers(np)
    fill!(local_grad_buffer, 0.0)

    Threads.@threads for j = 1:length(qdata)
        tid = Threads.threadid()

        iter_grad::Array{Float64} = iter_grad_buffer[:,:,tid]
        bv::BitVector = bv_buffer[:, tid]
        running_grad::Array{Float64} = running_grad_buffer[:, :, tid]
        
        fill!(iter_grad, 0.0)
        fill!(bv, false)
        fill!(running_grad, 1.0)

        eCF1::Float64, eCF2::Float64 = compute_eCF_and_gradient_recur!(qdata[j].eqn, params, iter_grad, bv, α, running_grad)
        eCF3::Float64 = 1 - eCF1 - eCF2

        eCF1 = max(eCF1, 1e-9)
        eCF2 = max(eCF2, 1e-9)
        eCF3 = max(eCF3, 1e-9)

        inv_eCF1::Float64 = 1 / eCF1
        inv_eCF2::Float64 = 1 / eCF2
        inv_eCF3::Float64 = 1 / eCF3

        total_loss_incr::Float64 = 
            ((q[j, 1] > 0) ? q[j, 1] * log(eCF1 / q[j, 1]) : 0.0) +
            ((q[j, 2] > 0) ? q[j, 2] * log(eCF2 / q[j, 2]) : 0.0) +
            ((q[j, 3] > 0) ? q[j, 3] * log(eCF3 / q[j, 3]) : 0.0)
        Threads.atomic_add!(total_loss, total_loss_incr)

        local_grad_buffer[:, tid] .+= q[j, 1] .* iter_grad[:, 1] .* inv_eCF1 +
            q[j, 2] .* iter_grad[:, 2] .* inv_eCF2 +
            q[j, 3] .* iter_grad[:, 3] .* inv_eCF3
    end

    gradient_storage .= sum(local_grad_buffer, dims=2)[:,1]
    return total_loss[]

end


"""
Recursive helper function that does the actual computations for [`compute_eCFs_and_gradient!`](@ref).
Returns eCFs for ab|cd and ac|bd -- ad|bc is calculated from the others.
"""
@fastmath function compute_eCF_and_gradient_recur!(
    eqn::RecursiveCFEquation, params::Vector{Float64},
    gradient_storage::Matrix{Float64},
    params_seen::BitVector,
    α::Float64,
    running_gradient::Array{Float64})::Tuple{Float64, Float64}

    if eqn.division_H == -1

        exp_sum = eqn.coal_edges == EMPTY_INT_VEC ? 1 : exp(-sum(params[j] for j = 1:length(params_seen) if eqn.coal_mask[j]))

        # Gradient computation
        @inbounds @simd for param_idx = 1:length(params)
            if eqn.coal_mask[param_idx]
                # we need to take the derivative
                @inbounds @simd for k = 1:3
                    gradient_storage[param_idx, k] += running_gradient[param_idx, k] * ((eqn.which_coal == k) ? 2/3*exp_sum : -1/3*exp_sum)
                end
            elseif params_seen[param_idx]
                # don't need to take the derivative, but we do need to
                # multiple the eCF value onto the gradient
                @inbounds @simd for k = 1:3
                    gradient_storage[param_idx, k] += running_gradient[param_idx, k] * ((eqn.which_coal == k) ? 1-2/3*exp_sum : 1/3*exp_sum)
                end
            end
        end

        # Return eCF contribution
        if eqn.which_coal == 1
            return 1-2/3*exp_sum, 1/3*exp_sum
        elseif eqn.which_coal == 2
            return 1/3*exp_sum, 1-2/3*exp_sum
        else
            return 1/3*exp_sum, 1/3*exp_sum
        end

    elseif length(eqn.divisions) == 4

        eqn_eCF1::Float64 = 0.0
        eqn_eCF2::Float64 = 0.0
        before = params_seen

        early_coal_exp_sum::Float64 = exp(-sum(params[j] for j = 1:length(params_seen) if eqn.coal_mask[j]))
        if eqn.can_coalesce_here
            # eCF contribution
            if eqn.which_coal == 1
                eqn_eCF1 += 1 - early_coal_exp_sum
            elseif eqn.which_coal == 2
                eqn_eCF2 += 1 - early_coal_exp_sum
            end

            # gradient contribution
            @inbounds @simd for param_idx = 1:length(params)
                if eqn.coal_mask[param_idx]
                    # calculate the derivative
                    !params_seen[param_idx] || error("Already seen this param??")
                    gradient_storage[param_idx, eqn.which_coal] += running_gradient[param_idx, eqn.which_coal] .* early_coal_exp_sum
                    params_seen[param_idx] = true
                elseif params_seen[param_idx]
                    gradient_storage[param_idx, eqn.which_coal] += running_gradient[param_idx, eqn.which_coal] .* (1 - early_coal_exp_sum)
                end
            end
        end

        running_gradient .*= early_coal_exp_sum
        !params_seen[eqn.division_H] || error("Already seen this param??")
        params_seen[eqn.division_H] = true
        γ::Float64 = params[eqn.division_H]
        @inbounds @simd for division_idx = 1:4
            split_grad::Float64 = quad_split_probability_gradient(division_idx, γ, α)
            split_prob::Float64 = quad_split_probability(division_idx, γ, α)
            prev_gamma_grad::Vector{Float64} = running_gradient[eqn.division_H, :]   # store here in case `split_grad` is 0.0

            # apply running gradient changes
            @inbounds @simd for param_idx = 1:length(params)
                if param_idx == eqn.division_H
                    running_gradient[eqn.division_H, :] .*= split_grad
                else
                    running_gradient[param_idx, :] .*= split_prob
                end
            end
            @inbounds @simd for e in eqn.coal_edges
                running_gradient[e, :] .*= -1.
            end

            recur_probs::Tuple{Float64, Float64} = compute_eCF_and_gradient_recur!(eqn.divisions[division_idx], params, gradient_storage, params_seen, α, running_gradient)

            # revert running gradient changes so that the next iteration is unbothered by them
            @inbounds @simd for e in eqn.coal_edges
                running_gradient[e, :] .*= -1.
            end
            @inbounds @simd for param_idx = 1:length(params)
                if param_idx == eqn.division_H
                    running_gradient[eqn.division_H, :] .= prev_gamma_grad
                else
                    running_gradient[param_idx, :] ./= split_prob
                end
            end

            recur_probs = early_coal_exp_sum * split_prob .* recur_probs

            eqn_eCF1 += recur_probs[1]
            eqn_eCF2 += recur_probs[2]
        end

        # revert running_gradient and params_seen changes
        running_gradient ./= early_coal_exp_sum
        @inbounds @simd for e in eqn.coal_edges
            params_seen[e] = false
        end
        params_seen[eqn.division_H] = false
        if !all(params_seen .== before)
            error("params_seen changed")
        end

        return eqn_eCF1, eqn_eCF2

    else

        !eqn.can_coalesce_here || error("Can coalesce w/ 1 taxa splitting at a single hybrid??")

        γ = params[eqn.division_H]
        params_seen[eqn.division_H] = true

        # apply running gradient changes
        @inbounds @simd for param_idx = 1:length(params)
            if param_idx != eqn.division_H
                running_gradient[param_idx, :] .*= γ
            end
        end

        eqn_eCF1, eqn_eCF2 = γ .* compute_eCF_and_gradient_recur!(eqn.divisions[1], params, gradient_storage, params_seen, α, running_gradient)

        # revert running gradient changes
        @inbounds @simd for param_idx = 1:length(params)
            if param_idx != eqn.division_H
                running_gradient[param_idx, :] ./= γ
            end
        end


        # apply running gradient changes
        @inbounds @simd for param_idx = 1:length(params)
            if param_idx != eqn.division_H
                running_gradient[param_idx, :] .*= (1-γ)
            else
                running_gradient[param_idx, :] .*= -1.
            end
        end

        secondary_probs = (1 - γ) .* compute_eCF_and_gradient_recur!(eqn.divisions[2], params, gradient_storage, params_seen, α, running_gradient)

        # revert running gradient changes
        @inbounds @simd for param_idx = 1:length(params)
            if param_idx != eqn.division_H
                running_gradient[param_idx, :] ./= (1-γ)
            else
                running_gradient[param_idx, :] .*= -1.
            end
        end

        params_seen[eqn.division_H] = false
        
        eqn_eCF1 += secondary_probs[1]
        eqn_eCF2 += secondary_probs[2]

        return eqn_eCF1, eqn_eCF2

    end
end


function quad_split_probability(type::Int, γ::Real, α::Real)::Float64
    if α == Inf
        # Strictly independent
        if type == 1
            return γ * γ
        elseif type == 2
            return (1 - γ) * (1 - γ)
        else
            return γ * (1 - γ)
        end
    elseif α == 0.0
        # Strictly dependent
        if type == 1
            return γ
        elseif type == 2
            return 1 - γ
        else
            return 0
        end
    else
        if type == 1
            # Same path, \gamma edge
            return γ * (1 / (α + 1) + α / (α + 1) * γ)
        elseif type == 2
            # Same path, 1-\gamma edge
            return (1 - γ) * (1 / (α + 1) + α / (α + 1) * (1 - γ))
        elseif type == 3 || type == 4
            # Different paths
            return γ * (α / (α + 1)) * (1 - γ)
        end
    end

    error("Found impossible type: $(type) (α = $(α))")
end


function quad_split_probability_gradient(type::Int, γ::Real, α::Real)::Float64
    if α == Inf
        # Strictly independent
        if type == 1
            return 2 * γ
        elseif type == 2
            return -2 * (1-γ)
        else
            return 1 - 2*γ
        end
    elseif α == 0.0
        # Strictly dependent
        if type == 1
            return 1
        elseif type == 2
            return -1
        else
            return 0.0
        end
    else
        # Correlated
        if type == 1
            return 1 / (α + 1) + 2 * γ * α / (α + 1)
        elseif type == 2
            return -1 / (α + 1) - 2 * (1 - γ) * α / (α + 1)
        elseif type == 3 || type == 4
            # return (α / (α + 1)) - 2 * γ * (α / (α + 1))
            return  (α / (α + 1)) * (1 - 2 * γ)
        end
    end

    error("Found impossible type: $(type) (α = $(α))")
end








