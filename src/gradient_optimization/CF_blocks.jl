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

@inline params_contained(block::DualGammaBlock) = [block.H]

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

@inline params_contained(block::EarlyCoalescenceBlock) = block.coal_edges


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

@inline params_contained(block::FailedEarlyCoalescenceBlock) = block.coal_edges


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

@inline params_contained(block::TwoTaxaHybridSplitBlock) = [block.H]


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

@inline params_contained(block::SimpleTreelikeBlock) = block.internal_edges


############ BLOCK PRODUCT ############

# This struct serves 1 purpose:
# To reduce the number of calls being made to `has_parameter` -- this
# function ALONE is taking up upwards of 20% of the `optimize_bls!` runtime!
struct BlockProduct
    blocks::AbstractArray{<:Block}
    relevant_parameters::AbstractArray{<:Int}   # array of param idxs that appear in this block product

    function BlockProduct(blocks::AbstractArray{<:Block})
        length(blocks) == 0 && new(blocks, [])
        new(blocks, union(reduce(vcat, params_contained(bl) for bl in blocks)))
    end
end


"""
Computes the value of a given block.
"""
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
Computes the derivative of a given block w.r.t. its parameters. Given that we are not taking
derivatives w.r.t. α, each block's derivative does not depend on the parameter it is taken
with respect to (so long as the block has the parameter of interest)
"""
function compute_block_deriv(block::Block, parameters::AbstractArray{<:Real}, eCF_type::Int, α::Real)
    if typeof(block) <: SimpleTreelikeBlock
        return compute_block_deriv(block, parameters, eCF_type)
    elseif typeof(block) <: TwoTaxaHybridSplitBlock
        return compute_block_deriv_wrt_γ(block, parameters, α)
    else
        return compute_block_deriv(block, parameters)
    end
end


function comp_all(bl, p)
    eCFs::Matrix{Float64} = zeros(size(bl))
    g = zeros(length(p))
    for j = 1:size(bl)[1]
        for k = 1:3
            compute_block_eCF_and_gradient!(bl[j,k], p, g, k, Inf)
        end
    end
    return nothing
end


"""
Computes the eCF and derivatives of a block at the same time to reduce redundant computations.
"""
function compute_block_eCF_and_gradient!(blocks::AbstractArray{<:AbstractArray{<:Block}}, parameters::AbstractArray{<:Real}, gradient_storage::AbstractArray{Float64}, eCF_type::Int, α::Real)
    
    block_values = fill(NaN, maximum(length(bv) for bv in blocks))
    block_derivs = fill(NaN, maximum(length(bv) for bv in blocks))
    
    eCF::Float64 = 0.0
    for block_prod in blocks
        block_derivs .= NaN

        # eCF computation
        vec_product::Float64 = 1.0
        for (block_idx, block) in enumerate(block_prod)
            block_values[block_idx] = compute_block_value(block, parameters, eCF_type, α)
            vec_product *= block_values[block_idx]
        end
        eCF += vec_product

        # gradient computation
        for param_idx = 1:length(parameters)
            blocks_have_param::Bool = false
            for bl in block_prod
                if has_parameter(bl, param_idx)
                    blocks_have_param = true
                    break
                end
            end
            blocks_have_param || continue
            
            param_deriv::Float64 = 1.0
            for (block_idx, block) in enumerate(block_prod)
                if has_parameter(block, param_idx)
                    if block_derivs[block_idx] === NaN
                        block_derivs[block_idx] = compute_block_deriv(block, parameters, eCF_type, α)
                    end
                    param_deriv *= block_derivs[block_idx]
                else
                    # These will not be NaN -- computed for the eCF calculation
                    param_deriv *= block_values[block_idx]
                end
            end
            gradient_storage[param_idx] += param_deriv
        end
    end
    return eCF
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

