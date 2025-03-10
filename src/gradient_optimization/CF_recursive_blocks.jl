

mutable struct RecursiveCFEquation
    can_coalesce_here::Bool
    coal_edges::AbstractArray{Edge}
    which_coal::Int     # 0 = NA, 1 = ab|cd, 2 = ac|bd, 3 = ad|bc
    division_H::Union{Nothing,Node}    # links to the hybrid node corresponding to γ
    divisions::AbstractArray{RecursiveCFEquation}
end


function get_blocks_from_recursive(eqn::RecursiveCFEquation)
    final_blocks_ab_cd = Vector{Vector{<:Block}}()
    final_blocks_ac_bd = Vector{Vector{<:Block}}()
    final_blocks_ad_bc = Vector{Vector{<:Block}}()
    
    get_blocks_from_recursive_recur!(eqn, final_blocks_ab_cd, Vector{Block}(), 1)
    get_blocks_from_recursive_recur!(eqn, final_blocks_ac_bd, Vector{Block}(), 2)
    get_blocks_from_recursive_recur!(eqn, final_blocks_ad_bc, Vector{Block}(), 3)

    return final_blocks_ab_cd, final_blocks_ac_bd, final_blocks_ad_bc
end

function get_blocks_from_recursive_recur!(eqn::RecursiveCFEquation, final_blocks::Vector{Vector{<:Block}}, contributing_blocks::Vector{<:Block}, which_quartet::Int)
    
    if eqn.division_H === nothing
        
        bl = SimpleTreelikeBlock(eqn.coal_edges, eqn.which_coal)
        push!(contributing_blocks, bl)
        push!(final_blocks, contributing_blocks)
    
    elseif length(eqn.divisions) == 4

        if eqn.can_coalesce_here && eqn.which_coal == which_quartet
            push!(final_blocks, [contributing_blocks; EarlyCoalescenceBlock(eqn.coal_edges)])
        end

        failed_early_block = FailedEarlyCoalescenceBlock(eqn.coal_edges)
        for division_idx = 1:4
            iter_bl = TwoTaxaHybridSplitBlock(eqn.division_H, division_idx)
            get_blocks_from_recursive_recur!(eqn.divisions[division_idx], final_blocks, [contributing_blocks; iter_bl; failed_early_block])
        end

    else    # length(eqn.divisions) == 2

        !eqn.can_coalesce_here || error("Can coalesce w/ 2 divisions??")
        for division_idx = 1:2
            iter_bl = DualGammaBlock(eqn.division_H, division_idx == 2)
            get_blocks_from_recursive_recur!(eqn.divisions[division_idx], final_blocks, [contributing_blocks; iter_bl], which_quartet)
        end

    end
        
end









