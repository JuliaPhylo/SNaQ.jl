using PhyloNetworks
const PN = PhyloNetworks;
const Node = PN.Node;
const Edge = PN.Edge;

abstract type Block end


############ DUAL GAMMA BLOCK ############

struct DualGammaBlock <: Block
    H::Int      # index of hybrid w/ the corresponding gamma value in `net.hybrid`
    major::Bool # did it take the major edge? i.e. equation is $\gamma$ if true, else $1-\gamma$
end

function compute_block_value(block::DualGammaBlock, params::AbstractArray{<:Real})::Float64
    return block.major ? 1 - params[block.H] : params[block.H]
end

function compute_block_deriv(block::DualGammaBlock, params::AbstractArray{<:Real})::Float64
    return block.major ? -1 : 1
end

has_parameter(block::DualGammaBlock, j::Int) = block.H == j


############ EARLY COALESCENCE BLOCK ############

struct EarlyCoalescenceBlock <: Block
    coal_edges::AbstractArray{Int}
end

function compute_block_value(block::EarlyCoalescenceBlock, params::AbstractArray{<:Real})::Float64
    return 1 - exp(-sum(params[e] for e in block.coal_edges))
end

function compute_block_deriv(block::EarlyCoalescenceBlock, params::AbstractArray{<:Real})::Float64
    return exp(-sum(params[e] for e in block.coal_edges))
end

has_parameter(block::EarlyCoalescenceBlock, j::Int) = j in block.coal_edges


############ FAILED EARLY COALESCENCE BLOCK ############

struct FailedEarlyCoalescenceBlock <: Block
    coal_edges::AbstractArray{Int}
end

function compute_block_value(block::FailedEarlyCoalescenceBlock, params::AbstractArray{<:Real})::Float64
    return exp(-sum(params[e] for e in block.coal_edges))
end

function compute_block_deriv(block::FailedEarlyCoalescenceBlock, params::AbstractArray{<:Real})::Float64
    return -exp(-sum(params[e] for e in block.coal_edges))
end

has_parameter(block::FailedEarlyCoalescenceBlock, j::Int) = j in block.coal_edges


############ TWO TAXA HYBRID SPLIT BLOCK ############

struct TwoTaxaHybridSplitBlock <: Block
    H::Int      # index of hybrid w/ the corresponding gamma value in `net.hybrid`
    type::Int   # 1: both taxa took the MINOR reticulation
                # 2: both taxa took the MAJOR reticulation
                # 3: lower taxa (according to `sort(.)`) took minor, other took major
                # 4: lower too major, other took minor
end

function compute_block_value(block::TwoTaxaHybridSplitBlock, params::AbstractArray{<:Real}, α::Real)::Float64
    γ = params[block.H]
    # @info "γ = $(γ), α = $(α)"
    
    if α == Inf
        # Strictly independent
        if block.type == 1
            return γ^2
        elseif block.type == 2
            return (1 - γ)^2
        else
            return γ * (1 - γ)
        end
    elseif α == 0.0
        # Strictly dependent
        if block.type == 1
            return γ
        elseif block.type == 2
            return 1 - γ
        else
            return 0
        end
    else
        if block.type == 1
            # @info "1: $(γ * (1 / (α + 1) + α / (α + 1) * γ))"
            return γ * (1 / (α + 1) + α / (α + 1) * γ)
        elseif block.type == 2
            # @info "2: $((1 - γ) * (1 / (α + 1) + α / (α + 1) * (1 - γ)))"
            return (1 - γ) * (1 / (α + 1) + α / (α + 1) * (1 - γ))
        elseif block.type == 3 || block.type == 4
            # @info "3: $(γ * (α / (α + 1)) * (1 - γ))"
            return γ * (α / (α + 1)) * (1 - γ)
        end
    end

    error("Found impossible block type: $(block.type) (α = $(α))")
end

function compute_block_deriv_wrt_γ(block::TwoTaxaHybridSplitBlock, params::AbstractArray{<:Real}, α::Real)::Float64
    γ = params[block.H]

    if α == Inf
        # Strictly independent
        if block.type == 1
            return 2 * γ
        elseif block.type == 2
            return -2 * (1-γ)
        else
            return 1 - 2*γ
        end
    elseif α == 0.0
        # Strictly dependent
        if block.type == 1
            return 1
        elseif block.type == 2
            return -1
        else
            return 0.0
        end
    else
        # Correlated
        if block.type == 1
            return 1 / (α + 1) + 2 * γ * α / (α + 1)
        elseif block.type == 2
            return -1 / (α + 1) - 2 * (1 - γ) * α / (α + 1)
        elseif block.type == 3 || block.type == 4
            # return (α / (α + 1)) - 2 * γ * (α / (α + 1))
            return  (α / (α + 1)) * (1 - 2 * γ)
        end
    end

    error("Found impossible block type: $(block.type) (α = $(α))")
end

has_parameter(block::TwoTaxaHybridSplitBlock, j::Int) = block.H == j


############ SIMPLE TREELIKE BLOCK ############

struct SimpleTreelikeBlock <: Block
    internal_edges::AbstractArray{Int}
    which_coal::Int     # 1: ab|cd
                        # 2: ac|bd
                        # 3: ad|bc
end

function compute_block_value(block::SimpleTreelikeBlock, params::AbstractArray{<:Real}, for_quartet::Int)::Float64
    if length(block.internal_edges) == 0 return 1/3 end
    exp_sum = exp(-sum(params[e] for e in block.internal_edges))
    return (block.which_coal == for_quartet) ? 1 - 2/3 * exp_sum : 1/3 * exp_sum
end

function compute_block_deriv(block::SimpleTreelikeBlock, params::AbstractArray{<:Real}, for_quartet::Int)::Float64
    if length(block.internal_edges) == 0 return 0.0 end
    exp_sum = exp(-sum(params[e] for e in block.internal_edges))
    return (block.which_coal == for_quartet) ? 2/3 * exp_sum : -1/3 * exp_sum
end

has_parameter(block::SimpleTreelikeBlock, j::Int) = j in block.internal_edges


"""
Computes the expected concordance factor that is defined by the list of multiplicative blocks in `blocks`.

`eCF_type` corresponds to which displayed quartet this eCF defines. I.e., `eCF_type=1` corresponds to ab|cd,
`2` corresponds to ac|bd, and `3` corresponds to ad|bc.
"""
function compute_eCF(blocks::AbstractArray{<:AbstractArray{<:Block}}, parameters::AbstractArray{<:Real}, eCF_type::Int, α::Real)::Float64
    eCF::Float64 = 0.0
    for block_vec in blocks
        vec_product::Float64 = 1.0
        for block in block_vec
            vec_product *= compute_block_value(block, parameters, eCF_type, α)
        end
        eCF += vec_product
    end
    return eCF
end

function compute_block_value(block::T, parameters::AbstractArray{<:Real}, eCF_type::Int, α::Real) where T <: Block
    if typeof(block) <: SimpleTreelikeBlock
        return compute_block_value(block, parameters, eCF_type)
    elseif typeof(block) <: TwoTaxaHybridSplitBlock
        return compute_block_value(block, parameters, α)
    else
        return compute_block_value(block, parameters)
    end
end


"""
Computes the derivative of the expected CF defined by `blocks` with respect to all parameters.
"""
function compute_block_derivs(blocks::AbstractArray{<:AbstractArray{<:Block}}, parameters::AbstractArray{<:Real}, eCF_type::Int, α::Real)
    
    block_values = fill(NaN, maximum(length(bv) for bv in blocks))
    block_derivs = fill(NaN, maximum(length(bv) for bv in blocks))

    derivs = zeros(length(parameters))
    for block_vec in blocks
        block_values .= NaN
        block_derivs .= NaN

        for param_idx = 1:length(parameters)
            blocks_have_param::Bool = false
            for bl in block_vec
                if has_parameter(bl, param_idx)
                    blocks_have_param = true
                    break
                end
            end
            blocks_have_param || continue
            
            param_deriv::Float64 = 1.0
            has_param::Bool = false
            for (block_idx, block) in enumerate(block_vec)
                if has_parameter(block, param_idx)
                    if block_derivs[block_idx] === NaN
                        block_derivs[block_idx] = compute_block_deriv(block, parameters, eCF_type, α)
                    end
                    param_deriv *= block_derivs[block_idx]
                else
                    if block_values[block_idx] === NaN
                        block_values[block_idx] = compute_block_value(block, parameters, eCF_type, α)
                    end
                    param_deriv *= block_values[block_idx]
                end
            end
            derivs[param_idx] += has_param ? param_deriv : 0.0
        end
    end
    return derivs
    
end

function compute_block_deriv(block::Block, parameters::AbstractArray{<:Real}, eCF_type::Int, α::Real)
    if typeof(block) <: SimpleTreelikeBlock
        return compute_block_deriv(block, parameters, eCF_type)
    elseif typeof(block) <: TwoTaxaHybridSplitBlock
        return compute_block_deriv_wrt_γ(block, parameters, α)
    else
        return compute_block_deriv(block, parameters)
    end
end




function compute_eCF(net::HybridNetwork, α::Real=Inf)

    param_map = Dict{Int, Int}()
    params = Array{Float64}(undef, length(net.hybrid) + length(net.edge))
    max_ID = net.numedges

    for (j, obj) in enumerate(vcat(net.hybrid, net.edge))
        obj.number = max_ID + j
        param_map[obj.number] = j
        params[j] = (typeof(obj) <: Node) ? getparentedgeminor(obj).gamma : obj.length
    end

    recur_eqns = get_reticulate_4taxa_quartet_equations(net, tipLabels(net), param_map)
    blocks1, blocks2, blocks3 = get_blocks_from_recursive(recur_eqns)

    eCFs = round.((
        compute_eCF(blocks1, params, 1, α),
        compute_eCF(blocks2, params, 2, α),
        compute_eCF(blocks3, params, 3, α)
    ), digits=5)
    @info eCFs

    return blocks1, blocks2, blocks3

end


function test_comp_eCF()
    n = readnewick("((a,b),(#H1,((c,d))#H1));")
    for E in n.edge
        if E.hybrid
            E.gamma = 0.5
        end
        E.length = 1.0
    end
    compute_eCF(n)  # (0.9773, 0.01135, 0.01135)
end


function test_comp_eCF2()
    n = readnewick("(((a,b),#H1),((c)#H1,d));")
    for E in n.edge
        if E.hybrid
            E.gamma = 0.5
        end
        E.length = 1.0
    end
    @info writenewick(n)
    compute_eCF(n)  # (0.86078, 0.06961, 0.06961)
end


function test_comp_eCF3()
    Random.seed!(42)
    n = readnewick("((a,(b,#H2)),(#H1,(((c)#H2,d))#H1));")
    for E in n.edge
        if E.hybrid && E.ismajor
            E.gamma = rand()
            E.gamma = max(E.gamma, 1 - E.gamma)
            getparentedgeminor(getchild(E)).gamma = 1 - E.gamma
        end
        E.length = rand()
    end
    @info writenewick(n)
    compute_eCF(n)  # (0.64646, 0.09303, 0.26051) - VERIFIED BY HAND
end


using Random
# seed=1: 
# seed=2: 
function test_comp_eCF_super_hard(; seed=abs(rand(Int)))
    @info seed
    Random.seed!(seed)
    n = readnewick("(((((a,#H5),((b,#H6))#H5))#H3,#H1),((#H2,((#H3,(((c)#H6,#H4),(d)#H4)))#H2))#H1);")
    for E in n.edge
        if E.hybrid && E.ismajor
            E.gamma = round(rand(), digits=3)
            E.gamma = max(E.gamma, 1 - E.gamma)
            getparentedgeminor(getchild(E)).gamma = 1 - E.gamma
        end
        E.length = 1.3 * round(rand(), digits=3)
    end
    @info writenewick(n)
    compute_eCF(n)
end

