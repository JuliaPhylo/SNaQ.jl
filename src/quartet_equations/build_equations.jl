# Building the equations for a whole set of quartets.

"""
Gathers a vector of `QuartetData` objects that define the expected
quartet concordance factors of `net`.
"""
findquartetequations(net::HybridNetwork)::Tuple{Vector{QuartetData},Dict,Vector{Float64},IdxObjMap,Vector{String}} =
    findquartetequations(net, 1:nchoose4taxalength(net))


 
"""
Deprecated - included for backwards compatibility in niche cases.
"""
find_quartet_equations(net::HybridNetwork) = findquartetequations(net)


"""
Gathers a vector of `QuartetData` objects that define the expected
quartet concordance factors of `net`. `q_idxs` is a `Vector{Int}` that
must be of length exactly (`net.numtaxa` choose 4). Each index of
`q_idxs` corresponds to a quartet whose equation will be computed.
"""
function findquartetequations(net::HybridNetwork, sampled_quartets::AbstractVector{Int})::Tuple{Vector{QuartetData},Dict,Vector{Float64},IdxObjMap,Vector{String}}
    fixnegativeedges!(net)
    all(e -> !e.hybrid || 1 >= e.gamma >= 0, net.edge) || error("net has gammas that are not in [0, 1]")
    all(h -> getparentedge(h).gamma + getparentedgeminor(h).gamma ≈ 1, net.hybrid) || error("net has hybrid with gammas that do not sum to 1")

    return findquartetequations!(net, sampled_quartets, Array{QuartetData}(undef, length(sampled_quartets)))
end


"""
See [`findquartetequations`](@ref)
"""
function findquartetequations!(net::HybridNetwork, sampled_quartets::AbstractVector{Int}, N_eqns::Vector{QuartetData})::Tuple{Vector{QuartetData},Dict,Vector{Float64},IdxObjMap,Vector{String}}
    # Relevant data to be returned
    t = sort(tiplabels(net))
    narg, param_map, idx_obj_map, params, _ = gatheroptimizationinfo(net)
    ntax::Int = length(t)
    nq::Int = length(sampled_quartets)
    ctx = treequartetcontext(net)   # built once for the whole batch, not per quartet
    scratch = pathscratchpool()     # ditto: one set of path buffers per thread

    if useparallelloop(nq)
        Threads.@threads for q_idx = 1:nq
            iter_taxa::AbstractVector{String} = t[unrank4taxa(ntax, sampled_quartets[q_idx])]
            N_eqns[q_idx] = findquartetequations4taxa(ctx, iter_taxa, param_map, 0.0,
                                                     scratch[Threads.threadid()])
        end
    else
        for q_idx = 1:nq
            iter_taxa = t[unrank4taxa(ntax, sampled_quartets[q_idx])]
            N_eqns[q_idx] = findquartetequations4taxa(ctx, iter_taxa, param_map, 0.0, scratch[1])
        end
    end

    return N_eqns, param_map, params, idx_obj_map, t
end


"""
    rebuildquartetequations!(net, old_eqns, new_eqns, ρ=0.0)

Recomputes every quartet equation for `net` into `new_eqns`, reusing the 4 taxa each entry
of `old_eqns` already records rather than unranking them again. The quartets a
[`search`](@ref) run uses are fixed for the whole run, so `old_eqns[j]` and `new_eqns[j]`
are always the same quartet.
"""
function rebuildquartetequations!(net::HybridNetwork, old_eqns::Vector{QuartetData},
                                  new_eqns::Vector{QuartetData}, ρ::Float64=0.0)
    param_map = gatheroptimizationinfo(net)[2]
    ctx = treequartetcontext(net)
    scratch = pathscratchpool()
    nq = length(old_eqns)
    if useparallelloop(nq)
        Threads.@threads for j = 1:nq
            new_eqns[j] = findquartetequations4taxa(ctx, old_eqns[j].q_taxa, param_map, ρ,
                                                   scratch[Threads.threadid()])
        end
    else
        for j = 1:nq
            new_eqns[j] = findquartetequations4taxa(ctx, old_eqns[j].q_taxa, param_map, ρ, scratch[1])
        end
    end
    return new_eqns
end


"""
[`findquartetequations4taxa`](@ref) via a precomputed [`TreeQuartetContext`](@ref), which
avoids re-deriving per quartet what is fixed for the whole network.
"""
function findquartetequations4taxa(ctx::TreeQuartetContext, taxa::AbstractVector{String}, parameter_map::Dict{Int, Int}, ρ::Float64=0.0,
                                   scratch::QuartetPathScratch=QuartetPathScratch())::QuartetData
    qdat = trytreelikequartet(ctx, taxa, parameter_map, scratch)
    qdat !== nothing && return qdat
    return reticulatequartetequations(ctx, taxa, parameter_map, ρ)
end


"""
Finds the quartet equations for the quarnet in `net` containing the taxa in `taxa`. `taxa` must contain exactly 4
    names of tips that are contained in `net`. `parameter_map` maps edges and gamma parameters in `net` to
    optimization variable indicies.
"""
function findquartetequations4taxa(net::HybridNetwork, taxa::AbstractVector{String}, parameter_map::Dict{Int, Int}, ρ::Float64=0.0)::QuartetData
    # Let's see if the quartet is tree-like and easy first
    qdat = trytreelikequartet(net, taxa, parameter_map)
    qdat !== nothing && return qdat
    return reticulatequartetequations(treequartetcontext(net), taxa, parameter_map, ρ)
end



