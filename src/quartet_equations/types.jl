# The expected-CF equation tree and the per-quartet data attached to it.

"""
    RecursiveCFEquation

The key struct used in computing -log pseudo-likelihood and gradients
when optimizing branch lengths.
"""
mutable struct RecursiveCFEquation
    can_coalesce_here::Bool
    coal_edges::Vector{Int}
    which_coal::Int     # 0 = NA, 1 = ab|cd, 2 = ac|bd, 3 = ad|bc
    division_H::Int     # hybrid index in `net.hybrid`; -1 means not set
    divisions::Vector{RecursiveCFEquation}
    # Dense `coal_edges` membership mask, built on first use by `coalmask`. The loss and
    # gradient read the sparse `coal_edges` directly, so it usually stays `nothing`.
    coal_mask::Union{Nothing,BitVector}
    nparam::Int

    function RecursiveCFEquation(can_coal::Bool, coal_Es::Vector{Int}, which_c::Int, dH::Int, d::Vector{RecursiveCFEquation}, nparam::Int)
        new(can_coal, coal_Es, which_c, dH, d, nothing, nparam)
    end
end


const EMPTY_EQN_VEC::Vector{RecursiveCFEquation} = Vector{RecursiveCFEquation}([]);


const EMPTY_INT_VEC::Vector{Int} = Vector{Int}([]);


"""
    coalmask(eqn::RecursiveCFEquation)::BitVector

`eqn`'s dense `coal_edges` membership mask, built and memoized on first request.
"""
function coalmask(eqn::RecursiveCFEquation)::BitVector
    m = eqn.coal_mask
    m === nothing || return m
    mask = falses(eqn.nparam)
    @inbounds for pidx in eqn.coal_edges
        mask[pidx] = true
    end
    eqn.coal_mask = mask
    return mask
end


"""
    QuartetData

A struct that contains:
1. The initial `RecursiveCFEquation` struct from which the loss & gradient can be calculated
2. A list of "internal" parameters (stored as indexed from 1 to `k` where `k` is the total number of
    parameters optimized during branch length optimization in `fitnumericalparameters!`) that are relevant to
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
