# Mapping a network onto the vector of numbers NLopt optimizes.

"""
Sets the branch lengths and γ values of edges in `net` according
to the values provided in `X` and the index-to-object map
provided in `idx_obj_map`.
"""
function setX!(net::HybridNetwork, X::Vector{Float64}, idx_obj_map::IdxObjMap)::Nothing
    for j in eachindex(X)
        obj::Union{Node, Edge} = idx_obj_map[j]
        if typeof(obj) <: PN.Node
            E_major = getparentedge(obj)
            E_minor = getparentedgeminor(obj)
            E_major.gamma = 1-X[j]
            E_minor.gamma = X[j]
        else
            obj.length = X[j]
        end
    end
end


"""
    countoptimizationparams(net)::Int

Number of parameters [`gatheroptimizationinfo`](@ref) would optimize for `net`, counted
directly. `gatheroptimizationinfo(net)[1]` answers the same question but builds the full
parameter map, bounds and starting values to do it -- and renumbers the network as a side
effect -- which is wasted work where only the count is wanted (e.g. deciding whether a
proposed move changed the parameter count).
"""
function countoptimizationparams(net::HybridNetwork)::Int
    n = length(net.hybrid)
    for e in net.edge
        childnode = getchild(e)
        childnode.leaf && continue
        if childnode.hybrid
            childnodechildren = getchildren(childnode)
            length(childnodechildren) == 1 && childnodechildren[1].leaf && continue
        end
        n += 1
    end
    return n
end


"""
Helper function to gather necessary information about `net` to
perform optimization.
"""
function gatheroptimizationinfo(net::HybridNetwork, change_numbers::Bool=true)
    param_map = Dict{Int, Int}()
    idx_obj_map::IdxObjMap = IdxObjMap();
    uq_ID = net.numedges
    param_idx = 1

    if change_numbers
        for obj in vcat(net.hybrid, net.edge, net.node)
            obj.number = uq_ID
            uq_ID += 1
        end
    end

    order = sortperm([obj.number for obj in vcat(net.hybrid, net.edge)])
    for obj in vcat(net.hybrid, net.edge)[order]
        if typeof(obj) <: Edge
            if getchild(obj).leaf continue end

            childnode = getchild(obj)
            if childnode.hybrid
                childnodechildren = getchildren(childnode)
                if length(childnodechildren) == 1 && childnodechildren[1].leaf
                    continue
                end
            end
            # if getchild(obj).hybrid && getchild(getchild(obj)).leaf continue end
        end

        haskey(param_map, obj.number) && error("Duplicate object number #$(obj.number).")
        param_map[obj.number] = param_idx
        idx_obj_map[param_idx] = obj
        param_idx += 1
    end

    params = gatherparams(net, param_map)
    narg = length(param_map)
    LB = Array{Float64}(undef, narg)
    UB = Array{Float64}(undef, narg)
    init_steps = Array{Float64}(undef, narg)

    for j = 1:narg
        obj::Union{Node,Edge} = idx_obj_map[j]
        if typeof(obj) <: Node
            params[j] = getparentedgeminor(obj).gamma
            LB[j] = 0.0
            UB[j] = 1.0
            init_steps[j] = 0.1
        else
            params[j] = obj.length
            LB[j] = 0.0
            UB[j] = 25.0
            init_steps[j] = 1.0
        end
    end

    return narg, param_map, idx_obj_map, params, LB, UB, init_steps
end


"""
Helper function that takes a network `net` and its `param_map` (provided
by [`gatheroptimizationinfo`](@ref)) and gathers each of the associated
parameters.
"""
function gatherparams(net::HybridNetwork, param_map::Dict{Int, Int})::Array{Float64}
    params = zeros(length(param_map))
    for obj in vcat(net.hybrid, net.edge)
        if haskey(param_map, obj.number)
            params[param_map[obj.number]] = typeof(obj) <: Node ? getparentedgeminor(obj).gamma : obj.length
        end
    end
    return params
end


"""
    gatherparams(net)

Helper function that takes a network `net` and its `param_map` (provided
by [`gatheroptimizationinfo`](@ref)) and gathers each of the associated
parameters.
"""
function gatherparams(net::HybridNetwork)::Array{Float64}
    param_map = gatheroptimizationinfo(net, true)[2]
    return gatherparams(net, param_map)
end
