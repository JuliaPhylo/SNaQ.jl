# Expected concordance factors and their gradients.

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
Sum of the branch lengths this equation coalesces over.
"""
@inline function coalbranchsum(eqn::RecursiveCFEquation, params::Vector{Float64})::Float64
    s = 0.0
    @inbounds for p in eqn.coal_edges
        s += params[p]
    end
    return s
end


"""
Recursive helper function that does the actual computations for [`computelossandgradient!`](@ref).
Returns eCFs for ab|cd and ac|bd -- ad|bc is calculated from the others. Gradient
contributions are added into `gradient_storage`, which the caller zeroes.

`rg` is the [`RunningGradient`](@ref) scratch space, which must be reset before the call and
is left reset afterwards: every exception added and every flag set is undone on the way out.
"""
@fastmath function computeexpectedCFandgradientrecur!(
    eqn::RecursiveCFEquation, params::Vector{Float64},
    gradient_storage::AbstractMatrix{Float64},
    rg::RunningGradient,
    α::Float64)::Tuple{Float64, Float64}

    coal_edges = eqn.coal_edges
    wc = eqn.which_coal

    if eqn.division_H == -1

        exp_sum = exp(-coalbranchsum(eqn, params))

        # Only two kinds of parameter get a contribution: the branches coalesced over here,
        # which take the derivative, and the ones an ancestor already fixed, which pick up
        # this block's probability. Both sets are small, so neither loop sweeps `params`.
        # d(eCF_k)/d(branch), the same for every branch coalesced over ...
        d1 = wc == 1 ? 2/3*exp_sum : -1/3*exp_sum
        d2 = wc == 2 ? 2/3*exp_sum : -1/3*exp_sum
        d3 = wc == 3 ? 2/3*exp_sum : -1/3*exp_sum
        # ... and eCF_k itself, which is what an already-fixed parameter picks up instead.
        f1 = wc == 1 ? 1-2/3*exp_sum : 1/3*exp_sum
        f2 = wc == 2 ? 1-2/3*exp_sum : 1/3*exp_sum
        f3 = wc == 3 ? 1-2/3*exp_sum : 1/3*exp_sum
        @inbounds for p in coal_edges
            gradient_storage[p, 1] += runninggradient(rg, p, 1) * d1
            gradient_storage[p, 2] += runninggradient(rg, p, 2) * d2
            gradient_storage[p, 3] += runninggradient(rg, p, 3) * d3
        end
        @inbounds for s = 1:rg.nexceptions
            p = rg.exceptionparam[s]
            rg.params_seen[p] || continue
            iscoaledge(coal_edges, p) && continue
            gradient_storage[p, 1] += rg.exceptionvalue[1, s] * f1
            gradient_storage[p, 2] += rg.exceptionvalue[2, s] * f2
            gradient_storage[p, 3] += rg.exceptionvalue[3, s] * f3
        end

        # Return eCF contribution
        if wc == 1
            return 1-2/3*exp_sum, 1/3*exp_sum
        elseif wc == 2
            return 1/3*exp_sum, 1-2/3*exp_sum
        else
            return 1/3*exp_sum, 1/3*exp_sum
        end

    elseif length(eqn.divisions) == 4

        eqn_eCF1::Float64 = 0.0
        eqn_eCF2::Float64 = 0.0
        nex0 = rg.nexceptions

        early_coal_exp_sum::Float64 = exp(-coalbranchsum(eqn, params))
        early_coal_exp_sum = max(early_coal_exp_sum, 1e-9)
        if eqn.can_coalesce_here
            # eCF contribution
            if wc == 1
                eqn_eCF1 += 1 - early_coal_exp_sum
            elseif wc == 2
                eqn_eCF2 += 1 - early_coal_exp_sum
            end

            # gradient contribution
            @inbounds for p in coal_edges
                rg.params_seen[p] && error("Already seen this param??")
                gradient_storage[p, wc] += runninggradient(rg, p, wc) * early_coal_exp_sum
            end
            @inbounds for s = 1:rg.nexceptions
                p = rg.exceptionparam[s]
                rg.params_seen[p] || continue
                gradient_storage[p, wc] += rg.exceptionvalue[wc, s] * (1 - early_coal_exp_sum)
            end
            @inbounds for p in coal_edges
                rg.params_seen[p] = true
            end
        end

        # Below this node the branches coalesced over here, and this reticulation's γ, no
        # longer scale like everything else, so they become exceptions.
        @inbounds for p in coal_edges
            addexception!(rg, p)
        end
        scalerunninggradient!(rg, early_coal_exp_sum)
        rg.params_seen[eqn.division_H] && error("Already seen this param??")
        addexception!(rg, eqn.division_H)
        rg.params_seen[eqn.division_H] = true
        γ::Float64 = params[eqn.division_H]

        for division_idx = 1:4
            split_grad::Float64 = quadsplitprobabilitygradient(division_idx, γ, α)
            split_prob::Float64 = quadsplitprobability(division_idx, γ, α)
            top = saverunninggradient!(rg)

            scalerunninggradient!(rg, split_prob, eqn.division_H)
            scaleexception!(rg, eqn.division_H, split_grad)
            @inbounds for e in coal_edges
                scaleexception!(rg, e, -1.0)
            end

            recur1::Float64, recur2::Float64 =
                computeexpectedCFandgradientrecur!(eqn.divisions[division_idx], params, gradient_storage, rg, α)

            restorerunninggradient!(rg, top)

            f = early_coal_exp_sum * split_prob
            eqn_eCF1 += f * recur1
            eqn_eCF2 += f * recur2
        end

        # revert running gradient and seen-flag changes
        scalerunninggradient!(rg, 1.0 / early_coal_exp_sum)
        @inbounds for e in coal_edges
            rg.params_seen[e] = false
        end
        rg.params_seen[eqn.division_H] = false
        dropexceptionsto!(rg, nex0)

        return eqn_eCF1, eqn_eCF2

    else

        !eqn.can_coalesce_here || error("Can coalesce w/ 1 taxa splitting at a single hybrid??")

        nex0 = rg.nexceptions
        γ = params[eqn.division_H]
        γ = max(1e-9, γ)
        γ = min(1.0 - 1e-9, γ)
        addexception!(rg, eqn.division_H)
        rg.params_seen[eqn.division_H] = true

        scalerunninggradient!(rg, γ, eqn.division_H)
        r1::Float64, r2::Float64 =
            computeexpectedCFandgradientrecur!(eqn.divisions[1], params, gradient_storage, rg, α)
        eqn_eCF1 = γ * r1
        eqn_eCF2 = γ * r2
        scalerunninggradient!(rg, 1.0 / γ, eqn.division_H)

        scalerunninggradient!(rg, 1 - γ, eqn.division_H)
        scaleexception!(rg, eqn.division_H, -1.0)
        s1::Float64, s2::Float64 =
            computeexpectedCFandgradientrecur!(eqn.divisions[2], params, gradient_storage, rg, α)
        scalerunninggradient!(rg, 1.0 / (1 - γ), eqn.division_H)
        scaleexception!(rg, eqn.division_H, -1.0)

        rg.params_seen[eqn.division_H] = false
        dropexceptionsto!(rg, nex0)

        eqn_eCF1 += (1 - γ) * s1
        eqn_eCF2 += (1 - γ) * s2

        return eqn_eCF1, eqn_eCF2

    end
end


"""
Whether `p` is one of the (few) branches in `coal`.
"""
@inline function iscoaledge(coal::Vector{Int}, p::Int)::Bool
    @inbounds for e in coal
        e == p && return true
    end
    return false
end


function quadsplitprobability(type::Int64, γ::Float64, α::Float64)::Float64
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


function quadsplitprobabilitygradient(type::Int64, γ::Float64, α::Float64)::Float64
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


function computeexpectedCF(qdata::QuartetData, params::Vector{Float64}, ρ::Float64=0.0)
    α = rhotoalpha(ρ)
    return computeexpectedCFandgradientrecur!(qdata.eqn, params, zeros(length(params), 3),
                                              RunningGradient(length(params)), α)
end


"""
    computeexpectedCF4taxa(net, taxa, ρ=0.0)

Helper function (primarily for the `QuartetNetworkGoodnessFit.jl` package) used
to compute the expected CFs of the quartet consisting of `taxa` in `net`.
The optional `ρ` argument (default 0) is the inheritance correlation parameter
in [0, 1]; `ρ = 0` is independent inheritance, `ρ = 1` is completely dependent.
"""
function computeexpectedCF4taxa(net::HybridNetwork, taxa::AbstractVector{<:AbstractString}, ρ::Real=0.0)::Tuple{Float64,Float64,Float64}
    # Definitely slightly inefficient to do this for each quartet, but shouldn't be a big deal.
    param_map, params = gatheroptimizationinfo(net)[[2,4]]

    qdata = findquartetequations4taxa(net, taxa, param_map)
    eCF1, eCF2 = computeexpectedCF(qdata, params, ρ)

    return eCF1, eCF2, 1-eCF1-eCF2
end


"""
    computeexpectedCFmatrix(net, ρ=0.0)

Computes the expected concordance factors of `net` with the inheritance
correlation parameter `ρ` (default=`0.0`). The returned `Matrix{Float64}`
object is unlabelled. See also [`computeexpectedDataCF`](@ref) for a `DataCF` object
with corresponding taxa information.
"""
function computeexpectedCFmatrix(net::HybridNetwork, ρ::Real=0.0)::Matrix{Float64}
    eqns, _, params, _ = findquartetequations(net)
    eCFs = zeros(length(eqns), 3)
    for j = 1:size(eCFs)[1]
        eCFs[j, 1], eCFs[j, 2] = computeexpectedCF(eqns[j], params, ρ)
        eCFs[j, 3] = 1 - eCFs[j, 1] - eCFs[j, 2]
    end
    return eCFs
end


"""
Deprecated - included for backwards compatibility in niche cases.
"""
computeexpectedCFs(net::HybridNetwork, ρ::Real=0.0)::Matrix{Float64} =
    computeexpectedCFmatrix(net, ρ)


"""
    computeexpectedDataCF(net, ρ=0.0)

Creates a DataCF object containing the expected CFs for each quartet in `net`.
"""
function computeexpectedDataCF(net::HybridNetwork, ρ::Real=0.0)::DataCF
    eqns, _, params, _ = findquartetequations(net);
    d = DataCF()
    for j in eachindex(eqns)
        eCF1, eCF2 = computeexpectedCF(eqns[j], params, ρ)
        q = Quartet(j, eqns[j].q_taxa..., Vector{Float64}([eCF1, eCF2, 1.0 - eCF1 - eCF2]))
        q.expCF = q.obsCF
        push!(d.quartet, q)
    end
    d.numQuartets = length(eqns)
    d.numTrees = -2 # code for expected CFs
    return d
end
