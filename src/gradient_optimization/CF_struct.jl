using Graphs    # for finding what edges contribute to internal edges; could remove at some point
using Combinatorics # for n choose 4, should definitely be removed later


struct eCFContribution
    hyb_edges::Vector{Int}
    internal_edges::Vector{Int}
    from_displayed_major::Bool
end

function find_quartet_equations(net::HybridNetwork)
    all(e -> e.length >= 0.0, net.edge) || error("net has negative edges")
    all(e -> !e.hybrid || 1 >= e.gamma >= 0, net.edge) || error("net has gammas that are not in [0, 1]")
    all(h -> getparentedge(h).gamma + getparentedgeminor(h).gamma ≈ 1, net.hybrid) || error("net has hybrid with gammas that do not sum to 1")

    t = sort(tipLabels(net))
    t_combos = combinations(t, 4)
    #quartet_eqns = Array{QEqn}(undef, numq)
    quartet_eqns = Array{Tuple{Vector{eCFContribution},Vector{eCFContribution},Vector{eCFContribution}}}(undef, length(t_combos))
    edge_number_to_idx_map = Dict{Int, Int}(e.number => j for (j, e) in enumerate(net.edge))
    quartet_taxa = Array{Vector{String}}(undef, length(t_combos))

    # Probably need to attach the taxanumbers index to quartet_eqns
    ts = [1,2,3,4]
    for j = 1:length(t_combos)
        quartet_eqns[j] = find_quartet_equations_4taxa(net, t[ts], edge_number_to_idx_map)
        quartet_taxa[j] = t[ts]

        ind = findfirst(x -> x>1, diff(ts))
        if ind === nothing ind = 4; end
        ts[ind] += 1
        for j in 1:(ind-1)
            ts[j] = j
        end
    end

    return quartet_eqns, t, quartet_taxa
end


function find_quartet_equations_4taxa(net::HybridNetwork, taxa::Vector{<:AbstractString}, edge_number_to_idx_map::Dict{Int, Int})
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

    # these will be our return values
    eCF12_34 = Vector{eCFContribution}()
    eCF13_24 = Vector{eCFContribution}()
    eCF14_23 = Vector{eCFContribution}()

    # iterate over the 2^H displayed trees
    G, W = Graph(net; withweights=true, minoredgeweight=Inf)
    for idx in eachindex(W) W[idx] = (W[idx] == Inf) ? Inf : 1.0 end
    H = net.hybrid
    node_to_idx = Dict{PN.Node, Int}(node => j for (j, node) in enumerate(net.node))                            # these two dicts used later for
    edge_to_graph_idxs = Dict{PN.Edge,Tuple{Int,Int}}(e => (node_to_idx[e.node[1]], node_to_idx[e.node[2]]) for e in net.edge) # easier graph weight adjustment
    H_edge_idx = Dict{PN.Edge, Int}(edge => j for (j, edge) in enumerate(net.edge) if edge.hybrid)  # maps each hybrid edge to its index
    H_minor_allowed = falses(length(H)) # each index corresponds to a hybrid in H
                                        # if the index is true, we are traversing the displayed
                                        # network that uses that hybrid's MAJOR parent,
                                        # else the one that uses that hybrid's MINOR parent
    
    for iter = 1:(2^length(H))
        if iter != 1 incr_bitvec!(H_minor_allowed) end
        
        # Enable/disable hybrid edges according to the BitVector
        currently_enabled_edges = Array{PN.Edge}(undef, length(H))
        for H_idx in 1:length(H)
            major_node1_idx, major_node2_idx = edge_to_graph_idxs[getparentedge(H[H_idx])]
            minor_node1_idx, minor_node2_idx = edge_to_graph_idxs[getparentedgeminor(H[H_idx])]
            if H_minor_allowed[H_idx]
                W[major_node1_idx, major_node2_idx] = W[major_node2_idx, major_node1_idx] = Inf
                W[minor_node1_idx, minor_node2_idx] = W[minor_node2_idx, minor_node1_idx] = 1.0
                currently_enabled_edges[H_idx] = getparentedgeminor(H[H_idx])
            else
                W[major_node1_idx, major_node2_idx] = W[major_node2_idx, major_node1_idx] = 1.0
                W[minor_node1_idx, minor_node2_idx] = W[minor_node2_idx, minor_node1_idx] = Inf
                currently_enabled_edges[H_idx] = getparentedge(H[H_idx])
            end
        end

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

            push!(eCF12_34, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                true
            ))
            push!(eCF13_24, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                false
            ))
            push!(eCF14_23, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                false
            ))
        elseif length(intersect(path_13, path_24)) == 0
            internal_graph_edges = intersect(path_12, path_34)    # 13|24 is displayed, so these paths must cross ONLY on the displayed edge portion
            internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)

            push!(eCF12_34, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                false
            ))
            push!(eCF13_24, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                true
            ))
            push!(eCF14_23, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                false
            ))
        else
            internal_graph_edges = intersect(path_12, path_34)    # 14|23 is displayed, so these paths must cross ONLY on the displayed edge portion
            internal_net_edges = from_graph_to_net_edges(net, internal_graph_edges)

            push!(eCF12_34, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                false
            ))
            push!(eCF13_24, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                false
            ))
            push!(eCF14_23, eCFContribution(
                [edge_number_to_idx_map[E.number] for E in currently_enabled_edges],
                [edge_number_to_idx_map[E.number] for E in internal_net_edges],
                true
            ))
        end
    end
    
    return eCF12_34, eCF13_24, eCF14_23
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
