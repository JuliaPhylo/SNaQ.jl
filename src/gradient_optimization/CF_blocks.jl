using PhyloNetworks
const PN = PhyloNetworks;
const Node = PN.Node;
const Edge = PN.Edge;

abstract type Block end


############ DUAL GAMMA BLOCK ############

struct DualGammaBlock <: Block
    H::Node     # hybrid w/ the corresponding gamma value
    major::Bool # did it take the major edge? i.e. equation is $\gamma$ if true, else $1-\gamma$
end

function compute_block_value(block::DualGammaBlock)::Float64
    return block.major ? getparentedge(block.H).gamma : getparentedgeminor(block.H).gamma
end

function compute_block_deriv(block::DualGammaBlock)::Float64
    return block.major ? -1 : 1
end

has_parameter(block::DualGammaBlock, e::Edge) = false
has_parameter(block::DualGammaBlock, H::Node) = block.H == H


############ EARLY COALESCENCE BLOCK ############

struct EarlyCoalescenceBlock <: Block
    coal_edges::AbstractArray{Edge}
end

function compute_block_value(block::EarlyCoalescenceBlock)::Float64
    return 1 - exp(-sum(e.length for e in block.coal_edges))
end

function compute_block_deriv(block::EarlyCoalescenceBlock)::Float64
    return exp(-sum(e.length for e in block.coal_edges)) 
end

has_parameter(block::EarlyCoalescenceBlock, e::Edge) = e in block.coal_edges
has_parameter(block::EarlyCoalescenceBlock, H::Node) = false


############ FAILED EARLY COALESCENCE BLOCK ############

struct FailedEarlyCoalescenceBlock <: Block
    coal_edges::AbstractArray{Edge}
end

function compute_block_value(block::FailedEarlyCoalescenceBlock)::Float64
    return exp(-sum(e.length for e in block.coal_edges))
end

function compute_block_deriv(block::FailedEarlyCoalescenceBlock)::Float64
    return -exp(-sum(e.length for e in block.coal_edges))
end

has_parameter(block::FailedEarlyCoalescenceBlock, e::Edge) = e in block.coal_edges
has_parameter(block::FailedEarlyCoalescenceBlock, H::Node) = false


############ TWO TAXA HYBRID SPLIT BLOCK ############

struct TwoTaxaHybridSplitBlock <: Block
    H::Node     # hybrid w/ the corresponding gamma
    type::Int   # 1: both taxa took the MINOR reticulation
                # 2: both taxa took the MAJOR reticulation
                # 3: lower taxa (according to `sort(.)`) took minor, other took major
                # 4: lower too major, other took minor
end

function compute_block_value(block::TwoTaxaHybridSplitBlock, α::Real)::Float64
    γ = getparentedgeminor(block.H).gamma
    
    if block.type == 1
        return γ * (1 / (α + 1) + α / (α + 1) * γ)
    elseif block.type == 2
        return (1 - γ) * (1 / (α + 1) + α / (α + 1) * (1 - γ))
    elseif block.type == 3 || block.type == 4
        return γ * (α / (α + 1)) * (1 - γ)
    else
        error("Found TwoTaxaHybridSplitBlock block.type = $(block.type)")
    end
end

function compute_block_deriv_wrt_γ(block::TwoTaxaHybridSplitBlock, α::Real)::Float64
    γ = getparentedgeminor(block.H).gamma

    if block.type == 1
        return 1 / (α + 1) + 2 * γ * α / (α + 1)
    elseif block.type == 2
        return -1 / (α + 1) - 2 * (1 - γ) * α / (α + 1)
    elseif block.type == 3 || block.type == 4
        # return (α / (α + 1)) - 2 * γ * (α / (α + 1))
        return  (α / (α + 1)) * (1 - 2 * γ)
    else
        error("Found TwoTaxaHybridSplitBlock block.type = $(block.type)")
    end
end

has_parameter(block::TwoTaxaHybridSplitBlock, e::Edge) = false
has_parameter(block::TwoTaxaHybridSplitBlock, H::Node) = block.H == H


############ SIMPLE TREELIKE BLOCK ############

struct SimpleTreelikeBlock <: Block
    internal_edges::AbstractArray{Edge}
    which_coal::Int     # 1: ab|cd
                        # 2: ac|bd
                        # 3: ad|bc
end

function compute_block_value(block::SimpleTreelikeBlock, for_quartet::Int)::Float64
    exp_sum = exp(-sum(e.length for e in block.internal_edges))
    return (block.which_coal == for_quartet) ? 1 - 2/3 * exp_sum : 1/3 * exp_sum
end

function compute_block_deriv(block::SimpleTreelikeBlock, for_quartet::Int)::Float64
    exp_sum = exp(-sum(e.length for e in block.internal_edges))
    return (block.which_coal == for_quartet) ? 2/3 * exp_sum : -1/3 * exp_sum
end

has_parameter(block::SimpleTreelikeBlock, e::Edge) = e in block.internal_edges
has_parameter(block::SimpleTreelikeBlock, H::Node) = false


"""
Computes the expected concordance factor that is defined by the list of multiplicative blocks in `blocks`.

`eCF_type` corresponds to which displayed quartet this eCF defines. I.e., `eCF_type=1` corresponds to ab|cd,
`2` corresponds to ac|bd, and `3` corresponds to ad|bc.
"""
function compute_eCF(blocks::AbstractArray{<:AbstractArray{<:Block}}, eCF_type::Int, edges::AbstractArray{Edge}, edge_number_to_idx_map::Dict{Int, Int}, α::Real)::Float64
    eCF::Float64 = 0.0
    for block_vec in blocks
        vec_product::Float64 = 1.0
        for block in block_vec
            vec_product *= compute_block_value(block, eCF_type, α)
        end
        eCF += vec_product
    end
    return eCF
end

function compute_block_value(block::T, eCF_type::Int, α::Real) where T <: Block
    if typeof(block) <: SimpleTreelikeBlock
        return compute_block_value(block, eCF_type)
    elseif typeof(block) <: TwoTaxaHybridSplitBlock
        return compute_block_value(block, α)
    else
        return compute_block_value(block)
    end
end


"""
Computes the derivative of the expected CF defined by `blocks` with respect to (wrt) either (a) an edge length
if `wrt` is of type `Edge`, or (b) a gamma value if `wrt` is of type `Node`.
"""
function compute_eCF_derivative(blocks::AbstractArray{<:AbstractArray{<:Block}}, eCF_type::Int, wrt::Union{Node,Edge}, edges::AbstractArray{Edge}, edge_number_to_idx_map::Dict{Int, Int}, α::Real)
    @warn "This function should likely be written the way the old one was so that we skip redundant work."
    
    deriv::Float64 = 0.0
    for block_vec in blocks
        param_in_vec = false
        vec_deriv::Float64 = 1.0
        for block in block_vec
            if !has_parameter(block, wrt)
                vec_deriv *= compute_block_value(block, eCF_type, α)
            else
                param_in_vec = true
                vec_deriv += compute_eCF_derivative(block, eCF_type, α)
            end
        end

        # If the parameter we are taking the derivative w.r.t. is not here, don't add to the derivative.
        if param_in_vec
            deriv += vec_deriv
        end
    end
    return deriv
end

function compute_eCF_derivative(block::Block, eCF_type::Int, α::Real)
    if typeof(block) <: SimpleTreelikeBlock
        return compute_block_deriv(block, eCF_type)
    elseif typeof(block) <: TwoTaxaHybridSplitBlock
        return compute_block_deriv_wrt_γ(block, α)
    else
        return compute_block_deriv(block)
    end
end






