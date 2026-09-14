# Equations for quartets a reticulation does span.

"""
    inducedquartetnetwork(ctx, taxa)

The sub-network of `ctx.net` spanned by the 4 leaves named in `taxa`: every node that is one
of them or an ancestor of one (following both parents at a reticulation), and every edge
between two such nodes. Node and edge numbers are preserved, since `parameter_map` is keyed
by them.

This is what pruning the network down to `taxa` produces, since deleting a leaf keeps exactly
the nodes that still have a descendant leaf, but it copies only the nodes involved rather
than the whole network.
"""
function inducedquartetnetwork(ctx::TreeQuartetContext, taxa::AbstractVector{String})::HybridNetwork
    net = ctx.net
    # Ancestral closure of the 4 leaves. Everything below is proportional to that closure
    # rather than to `net`, which is what `ctx`'s leaf and position maps are for.
    keep = Set{Int}()
    stack = Node[]
    kept = Node[]
    for name in taxa
        l = get(ctx.leafmap, name, nothing)
        l === nothing && error("inducedquartetnetwork: taxon $name is not a leaf of the network.")
        l.number in keep && continue
        push!(keep, l.number); push!(stack, l); push!(kept, l)
    end
    while !isempty(stack)
        n = pop!(stack)
        for p in getparents(n)
            p.number in keep && continue
            push!(keep, p.number)
            push!(stack, p)
            push!(kept, p)
        end
    end

    # `net.node` / `net.edge` order, so the result is deterministic (and identical to what
    # scanning the whole network in order would have produced).
    nodepos = ctx.nodepos
    sort!(kept, by = n -> nodepos[n.number])

    node_map = Dict{Int,Node}()
    sizehint!(node_map, length(kept))
    nodec = Vector{Node}(undef, length(kept))
    for (j, node) in enumerate(kept)
        nc = Node(node.number, node.leaf, node.hybrid)
        nc.name = node.name
        nodec[j] = nc
        node_map[node.number] = nc
    end

    # Every edge of the sub-network is incident to a kept node, so only their edge lists
    # need looking at; an edge with both endpoints kept shows up twice, hence `seen`.
    edges = Edge[]
    seen = Set{Int}()
    for node in kept, e in node.edge
        e.number in seen && continue
        n1 = e.node[1]; n2 = e.node[2]
        (haskey(node_map, n1.number) && haskey(node_map, n2.number)) || continue
        push!(seen, e.number)
        push!(edges, e)
    end
    edgepos = ctx.edgepos
    sort!(edges, by = e -> edgepos[e.number])

    edgec = Vector{Edge}(undef, length(edges))
    for (j, e) in enumerate(edges)
        c1 = node_map[e.node[1].number]; c2 = node_map[e.node[2].number]
        ec = Edge(e.number, e.length, e.hybrid, e.gamma, [c1, c2])
        ec.ischild1 = e.ischild1
        ec.ismajor = e.ismajor
        ec.containroot = e.containroot
        edgec[j] = ec
        push!(c1.edge, ec)
        push!(c2.edge, ec)
    end

    leafc = Node[nc for nc in nodec if nc.leaf]
    hybc = Node[nc for nc in nodec if nc.hybrid]
    rootnumber = net.node[net.rooti].number
    rooti = findfirst(nc -> nc.number == rootnumber, nodec)
    rooti === nothing && error("inducedquartetnetwork: the network root is not an ancestor of the quartet.")

    sub = HybridNetwork()
    sub.numtaxa = length(leafc)
    sub.numnodes = length(nodec)
    sub.numedges = length(edgec)
    sub.node = nodec
    sub.edge = edgec
    sub.leaf = leafc
    sub.rooti = rooti
    sub.names = net.names
    sub.hybrid = hybc
    sub.numhybrids = length(hybc)
    sub.isrooted = net.isrooted
    return sub
end


inducedquartetnetwork(net::HybridNetwork, taxa::AbstractVector{String})::HybridNetwork =
    inducedquartetnetwork(treequartetcontext(net), taxa)


"""
Quartet equations for a quartet a reticulation spans, so the tree-like shortcut does not
apply: prune down to `taxa` and split on the reticulations that remain. By far the most
expensive way to build one quartet's equations, and what `hmax > 0` runtime is dominated by.
"""
function reticulatequartetequations(ctx::TreeQuartetContext, taxa::AbstractVector{String}, parameter_map::Dict{Int, Int}, ρ::Float64=0.0)::QuartetData
    net = inducedquartetnetwork(ctx, taxa)

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
        get4taxaquartetequations(net, taxa, parameter_map, ρ),
        [parameter_map[obj.number] for obj in vcat(net.edge, net.hybrid)
            if haskey(parameter_map, obj.number)],
        taxa
    )
end


"""
Recursively builds the quartet CF equations for the quarnet in
    `net` containing the taxa in `taxa`. `taxa` should contain
    exactly 4 strings. `parameter_map` is a `Dict` that maps
    edge and hybrid node numbers (e.g. `node.number`) to a unique
    index - used for optimization. `ρ` is the inheritance correlation
    parameter.
"""
function get4taxaquartetequations(net::HybridNetwork, taxa::AbstractVector{String}, parameter_map::Dict{Int, Int}, ρ::Float64=0.0)::RecursiveCFEquation

    # If no hybrids remain, this case is simple
    if net.numhybrids == 0
        qdat = trytreelikequartet(net, taxa, parameter_map)
        return qdat.eqn
    end

    # Special case: we may still have hybrids, but 3+ leaves share a parent
    for L in net.leaf[[1, 2]]
        p_leaf = getparent(L)
        p_leaf_ch = getchildren(p_leaf)
        if length(p_leaf_ch) > 2 && sum(ch.leaf for ch in p_leaf_ch) > 2
            return RecursiveCFEquation(
                true, [], 1, -1, EMPTY_EQN_VEC, length(parameter_map)
            )
        end
    end


    # If >= 1 hybrids, we need to keep recursing
    lowest_H = getlowesthybrid(net)
    n_below_H = nleavesbelowlowesthybrid(lowest_H)
    if n_below_H == 1 && getparent(getparentedge(lowest_H)) != getparent(getparentedgeminor(lowest_H))
        # If there is only 1 leaf below this retic && the minor and major parents are the same,
        # we can just default to the final `else` case here and delete `lowest_H` b/c it has
        # no effect on eCFs

        # Remove the minor edge and all of its references in this copy
        div_major = deepcopynetwork(net)
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
        div_minor = deepcopynetwork(net)
        div_minor_H = div_minor.hybrid[findfirst(h -> h.name == lowest_H.name && h.number == lowest_H.number, div_minor.hybrid)]
        E_minor = getparentedgeminor(div_minor_H)
        E_major = getparentedge(div_minor_H)
        major_parent = getparent(E_major)
        for node in E_major.node
            node.edge = [e for e in node.edge if e != E_major]
        end
        PN.deleteEdge!(div_minor, E_major; part=false)
        PN.removeHybrid!(div_minor, getchild(E_minor))   # only removes its references - does not delete the node
        if length(getchildren(major_parent)) == 0
            major_parent.leaf = true
            push!(div_minor.leaf, major_parent)
            PN.deleteleaf!(div_minor, major_parent; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            major_parent.leaf = false
        end
        E_minor.hybrid = false
        E_minor.ismajor = true
        getchild(E_minor).hybrid = false


        # dmajnew = writenewick(div_major, round=true)
        # @info "DIV_MAJOR AFTER: $(dmajnew)"
        # dminnew = writenewick(div_minor, round=true)
        # @info "DIV_MINOR AFTER: $(dminnew)"
        r1 = get4taxaquartetequations(div_minor, taxa, parameter_map)
        r2 = get4taxaquartetequations(div_major, taxa, parameter_map)
        return RecursiveCFEquation(
            false, EMPTY_INT_VEC, 0, parameter_map[lowest_H.number],
            [r1, r2], length(parameter_map)
        )

    elseif n_below_H == 2
        # @info "4 - Following hybrid $(lowest_H.name)"
        #error("Implemented 2/4 cases where there are 2 leaves below hybrid so far - need to implement remaining 2 cases.")

        # @info net
        int_edges = getinternaledgesbelowlowesthybrid(lowest_H)
        leaves_below_H = getleavesbelowlowesthybrid(lowest_H)
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
        div1::HybridNetwork = deepcopynetwork(net)
        div1_H = div1.hybrid[findfirst(div1_H -> div1_H.number == lowest_H.number && div1_H.name == lowest_H.name, div1.hybrid)]
        E_minor = getparentedgeminor(div1_H)
        E_major = getparentedge(div1_H)

        # 1. Add placeholders for the new versions of the taxa
        #    and remove the current versions
        for L in leaves_below_H
            div1_L = div1.leaf[findfirst(dl -> dl.name == L.name, div1.leaf)]
            l = PN.addleaf!(div1, getchild(E_minor), "__$(L.name)", 0.0)
            PN.deleteleaf!(div1, div1_L; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            l.name = L.name
        end

        # 2. Delete hybrid edge - PhyloNetworks does all the clean up for us!
        PN.deletehybridedge!(div1, E_major, true, false, false, true, true)

        ######## Both taxa take the major edge ########
        # Same steps as above but for major instead of minor
        div2::HybridNetwork = deepcopynetwork(net)
        div2_H = div2.hybrid[findfirst(div2_H -> div2_H.number == lowest_H.number && div2_H.name == lowest_H.name, div2.hybrid)]
        E_minor = getparentedgeminor(div2_H)
        E_major = getparentedge(div2_H)

        # 1. Add placeholders for the new versions of the taxa
        #    and remove the current versions
        for L in leaves_below_H
            div2_L = div2.leaf[findfirst(dl -> dl.name == L.name, div2.leaf)]
            l = PN.addleaf!(div2, getchild(E_major), "__$(L.name)", 0.0)
            PN.deleteleaf!(div2, div2_L; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
            l.name = L.name
        end

        # 2. Delete hybrid edge - PhyloNetworks does all the clean up for us!
        PN.deletehybridedge!(div2, E_minor, true, false, false, true, true)


        ######## Lower taxa takes minor, higher takes major ########
        div3::HybridNetwork = deepcopynetwork(net)
        div3_H = div3.hybrid[findfirst(div3_H -> div3_H.number == lowest_H.number && div3_H.name == lowest_H.name, div3.hybrid)]
        E_minor = getparentedgeminor(div3_H)
        E_major = getparentedge(div3_H)

        # 1. add new version of lower leaf under minor retic's parent,
        #    then immediately delete the original leaf - PhyloNetworks
        #    takes care of net cleanup for us
        new_leaf::Node = PN.addleaf!(div3, getparent(E_minor), "__$(leaf_names[1])", 0.0)
        div3_L = div3.leaf[findfirst(dl -> dl.name == leaf_names[1], div3.leaf)]
        PN.deleteleaf!(div3, div3_L; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[1]

        # 2. vice versa
        new_leaf = PN.addleaf!(div3, getparent(E_major), "__$(leaf_names[2])", 0.0)
        div3_L = div3.leaf[findfirst(dl -> dl.name == leaf_names[2], div3.leaf)]
        PN.deleteleaf!(div3, div3_L; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[2]

        ######## Lower taxa takes major, higher takes minor ########
        div4::HybridNetwork = deepcopynetwork(net)
        div4_H = div4.hybrid[findfirst(div4_H -> div4_H.number == lowest_H.number && div4_H.name == lowest_H.name, div4.hybrid)]
        E_minor = getparentedgeminor(div4_H)
        E_major = getparentedge(div4_H)

        # 1. (same as above but flipped)
        new_leaf = PN.addleaf!(div4, getparent(E_major), "__$(leaf_names[1])", 0.0)
        div4_L = div4.leaf[findfirst(dl -> dl.name == leaf_names[1], div4.leaf)]
        PN.deleteleaf!(div4, div4_L; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[1]

        # 2. (same as above but flipped)
        new_leaf = PN.addleaf!(div4, getparent(E_minor), "__$(leaf_names[2])", 0.0)
        div4_L = div4.leaf[findfirst(dl -> dl.name == leaf_names[2], div4.leaf)]
        PN.deleteleaf!(div4, div4_L; simplify=false, nofuse=true, multgammas=false, keeporiginalroot=true)
        new_leaf.name = leaf_names[2]
        

        # @info "div1: $(writenewick(div1, round=true))"
        # @info "div2: $(writenewick(div2, round=true))"
        # @info "div3: $(writenewick(div3, round=true))"
        # @info "div4: $(writenewick(div4, round=true))"
        which_quartet = leaf_names[1] == taxa[1] ? (
            leaf_names[2] == taxa[2] ? 1 :
            leaf_names[2] == taxa[3] ? 2 : 3
        ) :
        leaf_names[1] == taxa[2] ? (
            leaf_names[2] == taxa[3] ? 3 : 2
        ) : 1
        recurrences::Array{RecursiveCFEquation} = Array{RecursiveCFEquation}(undef, 4)
        recurrences[1] = get4taxaquartetequations(div1, taxa, parameter_map)
        recurrences[2] = get4taxaquartetequations(div2, taxa, parameter_map)
        recurrences[3] = get4taxaquartetequations(div3, taxa, parameter_map)
        recurrences[4] = get4taxaquartetequations(div4, taxa, parameter_map)
        # @info "$(parameter_map[lowest_H.number]) -> $([eqn.division_H for eqn in recurrences])"

        return RecursiveCFEquation(
            length(int_edges) > 0, [parameter_map[int_e.number] for int_e in int_edges],
            which_quartet, parameter_map[lowest_H.number], recurrences, length(parameter_map)
        )
    else    # n_below_H is 3 or 4
        # 3 or 4 leaves below this hybrid, so it has no effect on eCFs!
        PN.deletehybridedge!(net, getparentedgeminor(lowest_H), false, true, false, true, false)    # params taken from blob deleting code
        return get4taxaquartetequations(net, taxa, parameter_map)
    end

end


"""
Gets the "lowest" hybrid, i.e. one of potentially multiple hybrids that do not have any other hybrids in their descendants.
Function assumes that extraneous retics have already been removed (i.e. retics on external quartet branches).
"""
function getlowesthybrid(net::HybridNetwork)::Node
    if net.numhybrids == 1 return net.hybrid[1] end
    return getlowesthybridrecur(net.hybrid[1])
end


"""
Helper function for [`getlowesthybrid`](@ref) - recursively finds the "lowest" hybrid in a network, starting at `node` - a hybrid node. 
"""
function getlowesthybridrecur(node::Node)
    if node.leaf
        return nothing
    end

    children = getchildren(node)
    for child in children
        child_val = getlowesthybridrecur(child)
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
function nleavesbelowlowesthybrid(H::Node)
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
function getinternaledgesbelowlowesthybrid(H::Node)::Vector{Edge}
    internal_edges = Vector{Edge}()
    c = getchildren(H)
    while length(c) == 1
        push!(internal_edges, getparentedge(c[1]))
        c = getchildren(c[1])
    end
    return internal_edges
end


"""
Helper function - gets the set of leaves below the hybrid node `H`.
"""
function getleavesbelowlowesthybrid(H::Node)::Vector{Node}
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
