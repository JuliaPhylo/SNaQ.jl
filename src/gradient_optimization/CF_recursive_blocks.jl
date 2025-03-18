


"""
The key struct used in computing -log pseudo-likelihood and gradients
when optimizing branch lengths.
"""
mutable struct RecursiveCFEquation
    can_coalesce_here::Bool
    coal_edges::AbstractArray{Int}
    which_coal::Int     # 0 = NA, 1 = ab|cd, 2 = ac|bd, 3 = ad|bc
    division_H::Int     # hybrid index in `net.hybrid`; -1 means not set
    divisions::AbstractArray{RecursiveCFEquation}
end


"""
TODO: remove (3) - it's just here for easier debugging

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
    relevant_params::AbstractArray{Int}
    q_taxa::Vector{String}
end
contains_parameter(qdata::QuartetData, param_idxs::AbstractVector{Int})::Bool = any(
    obj_idx -> obj_idx in qdata.relevant_params, param_idxs
)


function compute_eCF(qdata::QuartetData, params::Vector{<:Real})
    return compute_eCF_and_gradient_recur!(qdata.eqn, params, zeros(length(params), 3), falses(length(params)), Inf)
end


"""
Computes the loss (-log pseudo-likelihood) of network `N` given observed quartet concordance
factor data `q` under Dirichlet parameter `α`.
"""
function compute_loss(N::HybridNetwork, q, α::Real=Inf)::Float64
    N = deepcopy(N)
    semidirect_network!(N)
    qdata, _, params, _, _ = find_quartet_equations(N)
    return compute_loss(qdata, params, q, α)
end
function compute_loss(qdata::Vector{QuartetData}, params::AbstractVector{<:Real}, q, α::Real=Inf)::Float64
    return compute_loss_and_gradient!(qdata, params, zeros(length(params)), q, α)
end


"""
Computes expected concordance factors and gradients by recursively passing through `qdata`.
"""
function compute_loss_and_gradient!(qdata::Vector{QuartetData}, params::AbstractArray{<:Real}, gradient_storage::AbstractArray{Float64}, q, α::Real=Inf)::Float64

    bv::BitVector = falses(length(params))
    iter_grad::Array{Float64} = zeros(length(params), 3)
    fill!(gradient_storage, 0.0)

    total_loss::Float64 = 0.0
    for j = 1:length(qdata)
        fill!(iter_grad, 0.0)
        bv .= false

        eCF1, eCF2 = compute_eCF_and_gradient_recur!(qdata[j].eqn, params, iter_grad, bv, α)
        eCF3 = 1 - eCF1 - eCF2

        eCF1 = max(eCF1, 1e-9)
        eCF2 = max(eCF2, 1e-9)
        eCF3 = max(eCF3, 1e-9)

        total_loss += (q[j].data[1] > 0) ? q[j].data[1] * log(eCF1 / q[j].data[1]) : 0.0
        total_loss += (q[j].data[2] > 0) ? q[j].data[2] * log(eCF2 / q[j].data[2]) : 0.0
        total_loss += (q[j].data[3] > 0) ? q[j].data[3] * log(eCF3 / q[j].data[3]) : 0.0

        gradient_storage .+= q[j].data[1] .* iter_grad[:, 1] ./ eCF1
        gradient_storage .+= q[j].data[2] .* iter_grad[:, 2] ./ eCF2
        gradient_storage .+= q[j].data[3] .* iter_grad[:, 3] ./ eCF3

        # @info iter_grad[1,:]
        # @info round.([eCF1, eCF2, eCF3], digits=4)
        # @info round.(q[j].data[1:3], digits=4)
        # @info gradient_storage    # looking for 0.03319137891190929 on first iter in [1]
        #                           # 0.3143518750118932 on second iter in [1]
    end
    return total_loss

end


"""
Recursive helper function that does the actual computations for [`compute_eCFs_and_gradient!`](@ref).
Returns eCFs for ab|cd and ac|bd -- ad|bc is calculated from the others.
"""
function compute_eCF_and_gradient_recur!(
    eqn::RecursiveCFEquation, params::AbstractArray{<:Real},
    gradient_storage::AbstractArray{Float64},
    params_seen::BitVector,
    α::Real,
    running_gradient::Array{Float64}=ones(length(params), 3))::Tuple{Float64, Float64}

    # println("------------------------------------------------------------------------")
    # @info eqn
    # @info gradient_storage

    if eqn.division_H == -1

        # Simple treelike block
        exp_sum = length(eqn.coal_edges) == 0 ? 1 : exp(-sum(params[e] for e in eqn.coal_edges))

        # Gradient computation
        for param_idx = 1:length(params)
            if param_idx in eqn.coal_edges
                # we need to take the derivative
                for k = 1:3
                    gradient_storage[param_idx, k] += running_gradient[param_idx, k] * ((eqn.which_coal == k) ? 2/3*exp_sum : -1/3*exp_sum)
                end
            elseif params_seen[param_idx]
                # don't need to take the derivative, but we do need to
                # multiple the eCF value onto the gradient
                for k = 1:3
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
        before = deepcopy(params_seen)

        early_coal_exp_sum::Float64 = exp(-sum(params[e] for e in eqn.coal_edges))
        if eqn.can_coalesce_here
            # eCF contribution
            if eqn.which_coal == 1
                eqn_eCF1 += 1 - early_coal_exp_sum
            elseif eqn.which_coal == 2
                eqn_eCF2 += 1 - early_coal_exp_sum
            end

            # gradient contribution
            for param_idx = 1:length(params)
                if param_idx in eqn.coal_edges
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
        for division_idx = 1:4
            split_grad = quad_split_probability_gradient(division_idx, γ, α)
            split_prob = quad_split_probability(division_idx, γ, α)
            
            # apply running gradient changes
            for param_idx = 1:length(params)
                if param_idx == eqn.division_H
                    running_gradient[eqn.division_H, :] .*= split_grad
                else
                    running_gradient[param_idx, :] .*= split_prob
                end
            end
            for e in eqn.coal_edges
                running_gradient[e, :] .*= -1.
            end

            recur_probs = compute_eCF_and_gradient_recur!(eqn.divisions[division_idx], params, gradient_storage, params_seen, α, running_gradient)

            # revert running gradient changes so that the next iteration is unbothered by them
            for e in eqn.coal_edges
                running_gradient[e, :] .*= -1.
            end
            for param_idx = 1:length(params)
                if param_idx == eqn.division_H
                    running_gradient[eqn.division_H, :] ./= split_grad
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
        for e in eqn.coal_edges
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
        for param_idx = 1:length(params)
            if param_idx != eqn.division_H
                running_gradient[param_idx, :] .*= γ
            end
        end

        eqn_eCF1, eqn_eCF2 = γ .* compute_eCF_and_gradient_recur!(eqn.divisions[1], params, gradient_storage, params_seen, α, running_gradient)

        # revert running gradient changes
        for param_idx = 1:length(params)
            if param_idx != eqn.division_H
                running_gradient[param_idx, :] ./= γ
            end
        end


        # apply running gradient changes
        for param_idx = 1:length(params)
            if param_idx != eqn.division_H
                running_gradient[param_idx, :] .*= (1-γ)
            else
                running_gradient[param_idx, :] .*= -1.
            end
        end

        secondary_probs = (1 - γ) .* compute_eCF_and_gradient_recur!(eqn.divisions[2], params, gradient_storage, params_seen, α, running_gradient)

        # revert running gradient changes
        for param_idx = 1:length(params)
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
            return γ^2
        elseif type == 2
            return (1 - γ)^2
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
            # @info "1: $(γ * (1 / (α + 1) + α / (α + 1) * γ))"
            return γ * (1 / (α + 1) + α / (α + 1) * γ)
        elseif type == 2
            # @info "2: $((1 - γ) * (1 / (α + 1) + α / (α + 1) * (1 - γ)))"
            return (1 - γ) * (1 / (α + 1) + α / (α + 1) * (1 - γ))
        elseif type == 3 || type == 4
            # @info "3: $(γ * (α / (α + 1)) * (1 - γ))"
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








