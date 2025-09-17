# Includes code for updating `QuartetData` structs in-place based
# on what topological move was just conducted.
using PhyloNetworks


"""
Helper function to quickly determine whether or not an in-place update
exists for the move type `m`. Have this one function is cleaner and easier
to update when future in-place moves are implemented than by relying
on a bunch of `if` statements in various locations that might be missed.
"""
can_update_inplace(m::Symbol) = m in [:rNNI1, :rNNI2]


"""
Given that the move `move` was applied to network `Nprime` on parameters `params`, this
function updates `new_eqns` in-place with relevant changes to `old_eqns` where necessary.
"""
function update_quartet_equations!(
    old_eqns::Array{QuartetData},
    new_eqns::Array{QuartetData},
    Nprime::HybridNetwork,
    param_map::Dict{Int, Int},
    move::Symbol,
    params::Tuple,
    α::Real
)
    if move == :rNNI1
        apply_rNNI1_update!(Nprime, old_eqns, new_eqns, param_map, params[3], α)
    elseif move == :rNNI2
        apply_rNNI2_update!(Nprime, old_eqns, new_eqns, param_map, params..., α)
    else
        error("Only move that can be updated in place right now is rNNI1 (move = $(move))")
    end
end


function apply_rNNI1_update!(Nprime::HybridNetwork, old_qdata::AbstractVector{QuartetData}, new_qdata::AbstractVector{QuartetData}, param_map::Dict{Int, Int}, u::Node, α::Real=Inf)
    relevant_params = params_below_u_rNNI1(u, param_map)
    Threads.@threads for j in eachindex(old_qdata)
        if length(old_qdata[j].eqn.divisions) > 0 || recur_fxn_has_params(old_qdata[j].eqn, relevant_params)
            new_qdata[j] = find_quartet_equations_4taxa(Nprime, old_qdata[j].q_taxa, param_map, α)
        else
            new_qdata[j] = old_qdata[j]
        end
    end
end


function apply_rNNI2_update!(Nprime::HybridNetwork, old_qdata::AbstractVector{QuartetData}, new_qdata::AbstractVector{QuartetData}, param_map::Dict{Int, Int}, s::Node, t::Node, u::Node, v::Node, α::Real=Inf)
    relevant_params = [u.edge[findfirst(e -> t in e.node, u.edge)], s.edge[findfirst(e -> v in e.node, s.edge)]]
    relevant_params = [param_map[e.number] for e in relevant_params]
    Threads.@threads for j in eachindex(old_qdata)
        if contains_parameter(old_qdata[j], relevant_params)
            new_qdata[j] = find_quartet_equations_4taxa(Nprime, old_qdata[j].q_taxa, param_map, α)
        else
            new_qdata[j] = old_qdata[j]
        end
    end
end


"""
Gets the relevant objects (internal edges and γ's) underneath the node `u` where
`u` was used in an rNNI(1) move.
"""
function params_below_u_rNNI1(u::Node, param_map::Dict{Int,Int})::Vector{Int}
    objs_below = Vector{Int}([param_map[getparentedge(u).number]]);
    queue = Vector{Edge}([e for e in u.edge if getparent(e) == u]);
    n_iter::Int = 0

    while length(queue) > 0
        n_iter += 1
        n_iter < 1e7 || error("Stuck in infinite loop $(length(queue)), $(length(objs_below)).")

        next_edge = queue[length(queue)]
        deleteat!(queue, length(queue))
        
        if haskey(param_map, next_edge.number)
            if next_edge.number in objs_below
                continue
            else
                push!(objs_below, param_map[next_edge.number])
            end
        end
        child_node = getchild(next_edge)
        child_node.hybrid && !(child_node.number in objs_below) && push!(objs_below, param_map[child_node.number])
        for edge in child_node.edge
            edge != next_edge && getchild(edge) != child_node && push!(queue, edge)
        end

    end
    return objs_below
end


"""
Helpers function used to check whether any parameters in `param_idxs`
appear in the equation defined by `q`.
"""
function recur_fxn_has_params(q::RecursiveCFEquation, param_idxs::Vector{Int})::Bool
    any(e_idx -> e_idx in q.coal_edges || q.division_H === e_idx, param_idxs) && return true
    length(q.divisions) == 0 && return false
    return any(rec_q -> recur_fxn_has_params(rec_q, param_idxs), q.divisions)
end