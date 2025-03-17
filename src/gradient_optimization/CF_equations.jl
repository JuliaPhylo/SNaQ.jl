include("CF_struct.jl")
include("CF_blocks.jl")
include("CF_recursive_blocks.jl")
using PhyloNetworks, Graphs
const PN = PhyloNetworks;


function get_4taxa_quartet_equations(net::HybridNetwork, taxa::Vector{<:AbstractString}, parameter_map::Dict{Int, Int}, α::Real=Inf)::RecursiveCFEquation

    ########## REMOVE DEGREE-2 BLOBS ALONG EXTERNAL EDGES ##########
    # @info "BLOB BEFORE: $(writenewick(net, round=true))"
    # bcc = biconnectedcomponents(net, true) # true: ignore trivial blobs
    # entry = PN.biconnectedcomponent_entrynodes(net, bcc, true)
    # entryindex = indexin(entry, net.vec_node)
    # exitnodes = PN.biconnectedcomponent_exitnodes(net, bcc, false) # don't redo the preordering
    # bloborder = sortperm(entryindex) # pre-ordering for blobs in their own blob tree
    # function isexternal(ib) # is bcc[ib] of degree 2 and adjacent to an external edge?
    #     # yes if: 1 single exit adjacent to a leaf
    #     length(exitnodes[ib]) != 1 && return false
    #     ch = getchildren(exitnodes[ib][1])
    #     return length(ch) == 1 && ch[1].leaf
    # end
    # for ib in reverse(bloborder)
    #     isexternal(ib) || continue # keep bcc[ib] if not external of degree 2
    #     for he in bcc[ib]
    #         he.ismajor && continue
    #         # deletion of a hybrid can hide the deletion of another: check that he is still in net
    #         any(e -> e===he, net.edge) || continue
    #         # delete minor hybrid edge with options unroot=true: to make sure the
    #         # root remains of degree 3+, in case a degree-2 blob starts at the root
    #         # simplify=true: bc external blob
    #         PN.deletehybridedge!(net,he, false,true,false,true,false)
    #     end
    # end
    # @info "BLOB AFTER: $(writenewick(net, round=true))"
    ######################################################################

    # If no hybrids remain, this case is simple
    if net.numhybrids == 0
        quartet_type, int_edges = get_quartet_type_and_internal_edges(net, taxa, parameter_map)
        return RecursiveCFEquation(
            true, [parameter_map[int_e.number] for int_e in int_edges], quartet_type, -1,
            Vector{RecursiveCFEquation}([])
        )
    end

    # Special case: we may still have hybrids, but 3+ leaves share a parent
    for L in net.leaf[[1, 2]]
        p_leaf = getparent(L)
        p_leaf_ch = getchildren(p_leaf)
        if length(p_leaf_ch) > 2 && sum(ch.leaf for ch in p_leaf_ch) > 2
            return RecursiveCFEquation(
                true, [], 1, -1, Vector{RecursiveCFEquation}([])
            )
        end
    end


    # If >= 1 hybrids, we need to keep recursing
    lowest_H = get_lowest_hybrid(net)
    n_below_H = n_leaves_below_lowest_hybrid(lowest_H)
    if n_below_H == 1 && getparent(getparentedge(lowest_H)) != getparent(getparentedgeminor(lowest_H))
        # If there is only 1 leaf below this retic && the minor and major parents are the same,
        # we can just default to the final `else` case here and delete `lowest_H` b/c it has
        # no effect on eCFs


        # @info "2 - Following hybrid $(lowest_H.name)"

        # Remove the minor edge and all of its references in this copy
        div_major = deepcopy(net)
        div_major_H = div_major.hybrid[findfirst(h -> h.name == lowest_H.name && h.number == lowest_H.number, div_major.hybrid)]
        E_minor = getparentedgeminor(div_major_H)
        E_major = getparentedge(div_major_H)
        minor_parent = getparent(E_minor)
        for node in E_minor.node
            node.edge = [e for e in node.edge if e != E_minor]
        end
        PN.deleteEdge!(div_major, E_minor; part=false)
        PN.removeHybrid!(div_major, getchild(E_major))
        if length(getchildren(minor_parent)) == 0
            minor_parent.leaf = true
            push!(div_major.leaf, minor_parent)
            PN.deleteleaf!(div_major, minor_parent; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            minor_parent.leaf = false
        end
        E_major.hybrid = false
        E_major.ismajor = true
        getchild(E_major).hybrid = false

        # Remove the major edge and all of its references in this copy
        div_minor = deepcopy(net)
        div_minor_H = div_minor.hybrid[findfirst(h -> h.name == lowest_H.name && h.number == lowest_H.number, div_minor.hybrid)]
        E_minor = getparentedgeminor(div_minor_H)
        E_major = getparentedge(div_minor_H)
        major_parent = getparent(E_major)
        for node in E_major.node
            node.edge = [e for e in node.edge if e != E_major]
        end
        PN.deleteEdge!(div_minor, E_major; part=false)
        PN.removeHybrid!(div_minor, getchild(E_major))   # only removes its references - does not delete the node
        if length(getchildren(major_parent)) == 0
            major_parent.leaf = true
            push!(div_major.leaf, major_parent)
            PN.deleteleaf!(div_major, major_parent; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            major_parent.leaf = false
        end
        E_minor.hybrid = false
        E_minor.ismajor = true
        getchild(E_minor).hybrid = false

        # @info "DIV_MAJOR AFTER: $(writenewick(div_major, round=true))"
        # @info "DIV_MINOR: $(writenewick(div_minor, round=true))"
        r1 = get_4taxa_quartet_equations(div_minor, taxa, parameter_map)
        r2 = get_4taxa_quartet_equations(div_major, taxa, parameter_map)
        # @info "$(parameter_map[lowest_H.number]) -> $([eqn.division_H for eqn in [r1, r2]])"
        return RecursiveCFEquation(
            false, [], 0, parameter_map[lowest_H.number], [
                r1, r2
            ]
        )

    elseif n_below_H == 2
        # @info "4 - Following hybrid $(lowest_H.name)"
        #error("Implemented 2/4 cases where there are 2 leaves below hybrid so far - need to implement remaining 2 cases.")

        # @info net
        int_edges = get_internal_edges_below_lowest_hybrid(lowest_H)
        leaves_below_H = get_leaves_below_lowest_hybrid(lowest_H)
        leaf_names = sort([leaves_below_H[1].name, leaves_below_H[2].name])

        #### DEBUG STUFF##########################################################################
        # for (j, node) in enumerate(net.node)
        #     if node.name == ""
        #         node.name = "int$(j)"
        #     end
        # end
        # @info writenewick(net)
        # for E in int_edges
        #     @info "($(getparent(E).name), $(getchild(E).name))"
        # end
        ##########################################################################################

        ######## Both taxa take the minor edge ########
        div1 = deepcopy(net)
        div1_H = div1.hybrid[findfirst(div1_H -> div1_H.number == lowest_H.number && div1_H.name == lowest_H.name, div1.hybrid)]
        E_minor = getparentedgeminor(div1_H)
        E_major = getparentedge(div1_H)

        # 1. detach major edge, remove the hybrid status of the hybrid node
        for node in E_major.node
            node.edge = [e for e in node.edge if e != E_major]
        end
        PN.deleteEdge!(div1, E_major; part=false)
        PN.removeHybrid!(div1, getchild(E_minor))   # only removes its references - does not delete the node
        E_minor.hybrid = false
        E_minor.ismajor = true
        getchild(E_minor).hybrid = false

        # 2. for each of the 2 leaves, add a new leaf in the new location (under minor edge),
        #    and delete the original leaf immediately after. PhyloNetworks takes care of
        #    removing extraneous edges for us in the `PN.deleteleaf!` function.
        for L in leaves_below_H
            new_leaf = PN.addleaf!(div1, div1_H, "__$(L.name)", 0.0)
            PN.deleteleaf!(div1, L.name; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            new_leaf.name = L.name
        end

        ######## Both taxa take the major edge ########
        # Same steps as above but for major instead of minor
        div2 = deepcopy(net)
        div2_H = div2.hybrid[findfirst(div2_H -> div2_H.number == lowest_H.number && div2_H.name == lowest_H.name, div2.hybrid)]
        E_minor = getparentedgeminor(div2_H)
        E_major = getparentedge(div2_H)
        for node in E_minor.node
            node.edge = [e for e in node.edge if e != E_minor]
        end
        PN.deleteEdge!(div2, E_minor; part=false)
        PN.removeHybrid!(div2, getchild(E_major))   # only removes its references - does not delete the node
        E_major.hybrid = false
        getchild(E_major).hybrid = false

        for L in leaves_below_H
            new_leaf = PN.addleaf!(div2, div2_H, "__$(L.name)", 0.0)
            PN.deleteleaf!(div2, L.name; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            new_leaf.name = L.name
        end


        # Lower taxa takes minor, higher takes major
        div3 = deepcopy(net)
        div3_H = div3.hybrid[findfirst(div3_H -> div3_H.number == lowest_H.number && div3_H.name == lowest_H.name, div3.hybrid)]
        E_minor = getparentedgeminor(div3_H)
        E_major = getparentedge(div3_H)

        # 1. add new version of lower leaf under minor retic's parent,
        #    then immediately delete the original leaf - PhyloNetworks
        #    takes care of net cleanup for us
        new_leaf::Node = PN.addleaf!(div3, getparent(E_minor), "__$(leaf_names[1])", 0.0)
        PN.deleteleaf!(div3, leaf_names[1]; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[1]

        # 2. vice versa
        new_leaf = PN.addleaf!(div3, getparent(E_major), "__$(leaf_names[2])", 0.0)
        PN.deleteleaf!(div3, leaf_names[2]; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[2]

        # Lower taxa takes major, higher takes minor
        div4 = deepcopy(net)
        div4_H = div4.hybrid[findfirst(div4_H -> div4_H.number == lowest_H.number && div4_H.name == lowest_H.name, div4.hybrid)]
        E_minor = getparentedgeminor(div4_H)
        E_major = getparentedge(div4_H)

        # 1. (same as above but flipped)
        new_leaf = PN.addleaf!(div4, getparent(E_major), "__$(leaf_names[1])", 0.0)
        PN.deleteleaf!(div4, leaf_names[1]; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[1]

        # 2. (same as above but flipped)
        new_leaf = PN.addleaf!(div4, getparent(E_minor), "__$(leaf_names[2])", 0.0)
        PN.deleteleaf!(div4, leaf_names[2]; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[2]
        

        which_quartet = leaf_names[1] == taxa[1] ? (
            leaf_names[2] == taxa[2] ? 1 :
            leaf_names[2] == taxa[3] ? 2 : 3
        ) :
        leaf_names[1] == taxa[2] ? (
            leaf_names[2] == taxa[3] ? 3 : 2
        ) : 1
        recurrences = Array{RecursiveCFEquation}(undef, 4)
        recurrences[1] = get_4taxa_quartet_equations(div1, taxa, parameter_map)
        recurrences[2] = get_4taxa_quartet_equations(div2, taxa, parameter_map)
        recurrences[3] = get_4taxa_quartet_equations(div3, taxa, parameter_map)
        recurrences[4] = get_4taxa_quartet_equations(div4, taxa, parameter_map)
        # @info "$(parameter_map[lowest_H.number]) -> $([eqn.division_H for eqn in recurrences])"

        return RecursiveCFEquation(
            length(int_edges) > 0, [parameter_map[int_e.number] for int_e in int_edges],
            which_quartet, parameter_map[lowest_H.number], recurrences
        )
    else    # n_below_H is 3 or 4
        # 3 or 4 leaves below this hybrid, so it has no effect on eCFs!
        PN.deletehybridedge!(net, getparentedgeminor(lowest_H), false, true, false, true, false)    # params taken from blob deleting code
        return get_4taxa_quartet_equations(net, taxa, parameter_map)
    end

end


function get_quartet_type_and_internal_edges(net::HybridNetwork, taxa::Vector{<:AbstractString}, parameter_map::Dict{Int, Int})
    # iterate over the 2^H displayed trees
    G, W = Graph(net; withweights=true, minoredgeweight=Inf)
    for idx in eachindex(W) W[idx] = (W[idx] == Inf) ? Inf : 1.0 end
    node_to_idx = Dict{PN.Node, Int}(node => j for (j, node) in enumerate(net.node))                                            # these two dicts used later for
    edge_to_graph_idxs = Dict{PN.Edge,Tuple{Int,Int}}(e => (node_to_idx[e.node[1]], node_to_idx[e.node[2]]) for e in net.edge)  # easier graph weight adjustment
    
    # Find which edges form the internal edge of the quartet
    taxa1_node_idx = findfirst(n -> n.leaf && n.name == taxa[1], net.node)
    taxa2_node_idx = findfirst(n -> n.leaf && n.name == taxa[2], net.node)
    taxa3_node_idx = findfirst(n -> n.leaf && n.name == taxa[3], net.node)
    taxa4_node_idx = findfirst(n -> n.leaf && n.name == taxa[4], net.node)

    path_12 = a_star(G, taxa1_node_idx, taxa2_node_idx, W)
    path_34 = a_star(G, taxa3_node_idx, taxa4_node_idx, W)
    path_13 = a_star(G, taxa1_node_idx, taxa3_node_idx, W)
    path_24 = a_star(G, taxa2_node_idx, taxa4_node_idx, W)

    # edges here have direction, but we don't want them to, so we set (src,dst) = (minidx,maxidx) so the direction is effectively removed
    for path in [path_12, path_34, path_13, path_24]
        for (j, edge) in enumerate(path)
            path[j] = Graphs.SimpleEdge(min(edge.src, edge.dst), max(edge.src, edge.dst))
        end
    end

    if length(intersect(path_12, path_34)) == 0
        internal_graph_edges = intersect(path_13, path_24)    # 12|34 is displayed, so these paths must cross ONLY on the displayed edge portion
        internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)
        return 1, internal_net_edges
    elseif length(intersect(path_13, path_24)) == 0
        internal_graph_edges = intersect(path_12, path_34)    # 13|24 is displayed, so these paths must cross ONLY on the displayed edge portion
        internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)
        return 2, internal_net_edges
    else
        internal_graph_edges = intersect(path_12, path_34)    # 14|23 is displayed, so these paths must cross ONLY on the displayed edge portion
        internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)
        return 3, internal_net_edges
    end
end


"""
Gets the "lowest" hybrid, i.e. one of potentially multiple hybrids that do not have any other hybrids in their descendants.
Function assumes that extraneous retics have already been removed (i.e. retics on external quartet branches).
"""
function get_lowest_hybrid(net::HybridNetwork)::Node
    if net.numhybrids == 1 return net.hybrid[1] end
    return get_lowest_hybrid_recur(net.hybrid[1])
end

function get_lowest_hybrid_recur(node::Node)
    if node.leaf
        return nothing
    end

    children = getchildren(node)
    for child in children
        child_val = get_lowest_hybrid_recur(child)
        if child_val !== nothing return child_val end
    end

    if node.hybrid
        return node
    else
        return nothing
    end
end


"""
Gets the number of leaves in a quarnet below the lowest hybrid in the quarnet.
Assumes that reticulations on external edges are removed. HOWEVER there may
still be more than 2 leaves below a hybrid.
"""
function n_leaves_below_lowest_hybrid(H::Node)
    queue = getchildren(H)
    leaves_found::Int = 0
    while length(queue) > 0
        curr = queue[length(queue)]
        deleteat!(queue, length(queue))
        if curr.leaf
            leaves_found += 1
        else
            for c in getchildren(curr)
                push!(queue, c)
            end
        end
    end
    return leaves_found
end


"""
Assumes that there are 2 leaves below `H` in the quarnet.
"""
function get_internal_edges_below_lowest_hybrid(H::Node)
    internal_edges = Vector{Edge}()
    c = getchildren(H)
    while length(c) == 1
        push!(internal_edges, getparentedge(c[1]))
        c = getchildren(c[1])
    end
    return internal_edges
end


function get_leaves_below_lowest_hybrid(H::Node)
    queue = Vector{Node}([H])
    leaves = Vector{Node}([])

    while length(queue) > 0
        curr = queue[length(queue)]
        deleteat!(queue, length(queue))

        if curr.leaf push!(leaves, curr) end
        for c in getchildren(curr)
            push!(queue, c)
        end
    end

    return leaves
end


function find_quartet_equations(net::HybridNetwork, recur_eqns::Array{QuartetData}=Vector{QuartetData}([]))
    all(e -> e.length >= 0.0, net.edge) || error("net has negative edges")
    all(e -> !e.hybrid || 1 >= e.gamma >= 0, net.edge) || error("net has gammas that are not in [0, 1]")
    all(h -> getparentedge(h).gamma + getparentedgeminor(h).gamma ≈ 1, net.hybrid) || error("net has hybrid with gammas that do not sum to 1")

    t = sort(tipLabels(net))
    t_combos = combinations(t, 4)
    #quartet_eqns = Array{QEqn}(undef, numq)

    narg, param_map, idx_obj_map, params, _ = gather_optimization_info(net)

    if length(recur_eqns) == 0
        recur_eqns = Array{QuartetData}(undef, length(t_combos))
    end

    # Probably need to attach the taxanumbers index to quartet_eqns
    ts = [1,2,3,4]
    for j = 1:length(t_combos)
        # reduced_net = get_reduced_net(reduced_nets, Set{String}(t[ts]), t)
        # recur_eqns[j] = find_quartet_equations_4taxa(reduced_net, t[ts], param_map)
        recur_eqns[j] = find_quartet_equations_4taxa(net, t[ts], param_map)     # 625ms, 287MiB

        incr_taxa_idx!(ts)
    end

    return recur_eqns, param_map, params, idx_obj_map, t
end


"""
Helper function to increment the 4-taxa index within `find_quartet_equations`.
"""
function incr_taxa_idx!(ts::Vector{Int})
    ind = findfirst(x -> x>1, diff(ts))
    if ind === nothing ind = 4; end
    ts[ind] += 1
    for j in 1:(ind-1)
        ts[j] = j
    end
end


"""
Dynamic programming approach to iteratively reducing the initial hybrid network down to quarnets.
"""
function get_reduced_net(reduced_nets::Dict{Set{String}, HybridNetwork}, taxa::Set{String}, all_taxa::Vector{String})::HybridNetwork
    haskey(reduced_nets, taxa) && return reduced_nets[taxa]

    next_taxa::String = all_taxa[findfirst(t -> !(t in taxa), all_taxa)]
    push!(taxa, next_taxa)
    parent_net::HybridNetwork = get_reduced_net(reduced_nets, taxa, all_taxa)
    
    net::HybridNetwork = deepcopy(parent_net)
    PN.deleteleaf!(net, next_taxa, simplify=false, unroot=false, nofuse=true)
    delete!(taxa, next_taxa)
    net.numtaxa == length(taxa) || error("$(net.numtaxa) != $(length(taxa))")

    # find and delete degree-2 blobs along external edges
    bcc = biconnectedcomponents(net, true) # true: ignore trivial blobs
    entry = PN.biconnectedcomponent_entrynodes(net, bcc, true)
    entryindex = indexin(entry, net.vec_node)
    exitnodes = PN.biconnectedcomponent_exitnodes(net, bcc, false) # don't redo the preordering
    bloborder = sortperm(entryindex) # pre-ordering for blobs in their own blob tree
    function isexternal(ib) # is bcc[ib] of degree 2 and adjacent to an external edge?
        # yes if: 1 single exit adjacent to a leaf
        length(exitnodes[ib]) != 1 && return false
        ch = getchildren(exitnodes[ib][1])
        return length(ch) == 1 && ch[1].leaf
    end
    for ib in reverse(bloborder)
        isexternal(ib) || continue # keep bcc[ib] if not external of degree 2
        for he in bcc[ib]
            he.ismajor && continue
            # deletion of a hybrid can hide the deletion of another: check that he is still in net
            any(e -> e===he, net.edge) || continue
            # delete minor hybrid edge with options unroot=true: to make sure the
            # root remains of degree 3+, in case a degree-2 blob starts at the root
            # simplify=true: bc external blob
            PN.deletehybridedge!(net,he, false,true,false,true,false)
        end
    end
    reduced_nets[deepcopy(taxa)] = net

    return net
end


function find_quartet_equations_4taxa(net::HybridNetwork, taxa::Vector{<:AbstractString}, parameter_map::Dict{Int, Int}, α::Real=Inf)::QuartetData
    net = deepcopy(net)
    
    # remove all taxa other than those in `taxa`
    for t in sort(tipLabels(net))
        t in taxa && continue
        PhyloNetworks.deleteleaf!(net, t, simplify=false, unroot=false, nofuse=true)
    end

    # find and delete degree-2 blobs along external edges
    bcc = biconnectedcomponents(net, true) # true: ignore trivial blobs
    entry = PN.biconnectedcomponent_entrynodes(net, bcc, true)
    entryindex = indexin(entry, net.vec_node)
    exitnodes = PN.biconnectedcomponent_exitnodes(net, bcc, false) # don't redo the preordering
    bloborder = sortperm(entryindex) # pre-ordering for blobs in their own blob tree
    function isexternal(ib) # is bcc[ib] of degree 2 and adjacent to an external edge?
        # yes if: 1 single exit adjacent to a leaf
        length(exitnodes[ib]) != 1 && return false
        ch = getchildren(exitnodes[ib][1])
        return length(ch) == 1 && ch[1].leaf
    end
    for ib in reverse(bloborder)
        isexternal(ib) || continue # keep bcc[ib] if not external of degree 2
        for he in bcc[ib]
            he.ismajor && continue
            # deletion of a hybrid can hide the deletion of another: check that he is still in net
            any(e -> e===he, net.edge) || continue
            # delete minor hybrid edge with options unroot=true: to make sure the
            # root remains of degree 3+, in case a degree-2 blob starts at the root
            # simplify=true: bc external blob
            PN.deletehybridedge!(net,he, false,true,false,true,false)
        end
    end

    return QuartetData(
        get_4taxa_quartet_equations(net, taxa, parameter_map, α),
        [parameter_map[obj.number] for obj in vcat(net.edge, net.hybrid)
            if haskey(parameter_map, obj.number)],
        taxa
    )
end


function compute_logPL(N::HybridNetwork, obsCFs)::Float64
    eqns, _, params, _ = find_quartet_equations(N)
    total_loss = 0.0
    for j = 1:size(eqns)[1]
        for k = 1:3
            # Loss function
            total_loss += obsCFs[j].data[k] * log(compute_eCF(eqns[j, k], params, k, Inf) / obsCFs[j].data[k])
        end
    end
    return total_loss
end



function from_graph_to_net_edges(net::HybridNetwork, internal_graph_edges::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}})
    net_edges = Array{PN.Edge}(undef, length(internal_graph_edges))
    for (E_idx, E) in enumerate(internal_graph_edges)
        nodei = net.node[E.src]
        nodej = net.node[E.dst]
        edgeij = nodei.edge[findfirst(e -> nodej in e.node, nodei.edge)]
        net_edges[E_idx] = edgeij
    end
    return net_edges
end


"""
Function taken from InPhyNet.jl

Converts the tree/network `net` into a SimpleGraph to leverage already
implemented pathfinding algorithms.

# Arguments
- includeminoredges (default=true): if true, the entire network is translated to a graph.
      Otherwise, only tree-like edges (other than those in `alwaysinclude`) are retained.
- alwaysinclude (default=nothing): edges that should always be included in the graph,
      regardless of the value of `includeminoredges`
- withweights (default=false): return a set of weights corresponding to branch lengths as well
"""
function Graph(net::HybridNetwork; withweights::Bool=false, minoredgeweight::Float64=1.)
    graph = SimpleGraph(net.numnodes)
    weights = Matrix{Float64}(undef, net.numnodes, net.numnodes)
    weights .= Inf
    nodemap = Dict{Node, Int64}(node => idx for (idx, node) in enumerate(net.node))
    for edge in net.edge
        enode1 = edge.node[1]
        enode2 = edge.node[2]
        if haskey(nodemap, enode1) && haskey(nodemap, enode2)
            add_edge!(graph, nodemap[enode1], nodemap[enode2])
            if withweights
                weight = 1
                if edge.hybrid && !edge.ismajor
                    weight = minoredgeweight
                elseif edge.length == -1.
                    weight = 1
                end
                
                weights[nodemap[enode1], nodemap[enode2]] =
                    weights[nodemap[enode2], nodemap[enode1]] = weight
            end
        end
    end

    if withweights
        return graph, weights
    end
    return graph
end
