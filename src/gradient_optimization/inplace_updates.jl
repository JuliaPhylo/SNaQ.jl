# Includes code for updating `QuartetData` structs in-place based
# on what topological move was just conducted.
using PhyloNetworks


function apply_rNNI1_update!(Nprime::HybridNetwork, old_qdata::AbstractVector{QuartetData}, new_qdata::AbstractVector{QuartetData}, param_map::Dict{Int, Int}, u::Node, α::Real=Inf)
    relevant_params = params_below_u_rNNI1(u, param_map)
    for j = 1:length(old_qdata)
        if recur_fxn_has_params(old_qdata[j].eqn, relevant_params)
            new_qdata[j] = find_quartet_equations_4taxa(Nprime, old_qdata[j].q_taxa, param_map, α)
        else
            new_qdata[j] = old_qdata[j]
        end
    end
end


function apply_rNNI2_update!(Nprime::HybridNetwork, old_qdata::AbstractVector{QuartetData}, new_qdata::AbstractVector{QuartetData}, param_map::Dict{Int, Int}, s::Node, t::Node, u::Node, v::Node, α::Real=Inf)
    relevant_params = [u.edge[findfirst(e -> t in e.node, u.edge)], s.edge[findfirst(e -> v in e.node, s.edge)]]
    relevant_params = [param_map[e.number] for e in relevant_params]
    for j = 1:length(old_qdata)
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