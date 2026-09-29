# General network helpers shared by the inference code.

"""
    fixnegativeedges!(N::HybridNetwork)

Sets negative edge lengths to 0.0 if they do not lead to leaves. These are
most often hybrid edges that were not given edge lengths because their
edge lengths were unidentifiable.
"""
function fixnegativeedges!(N::HybridNetwork; shouldwarn::Bool=true)
    for E in N.edge
        getchild(E).leaf && continue
        if E.hybrid && E.length < 0.0
            shouldwarn && @warn "Found a hybrid edge with unspecified length; setting its length to 0.0 before proceeding. To avoid this, set the edge's length manually."
            shouldwarn = false
            E.length = 0.0
        end

        E.length < 0 && error("""
        Found an edge with unspecified length that was neither a hybrid edge nor a leaf edge. This is not allowed.
        This is likely due to the edge not having a specified length in the newick string used to read the network.
        If this is not the case, please submit a bug report at github.com/JuliaPhylo/SNaQ.jl/issues

        To identify such edges in your network, run the following (after changing `net` to your network's variable):
            badedges = getnegativeedges(net)
        """)
    end
end


"""
    getnegativeedges(N::HybridNetwork)

Helper function to get all the edges in the network `N` that have unspecified lengths.
"""
function getnegativeedges(N::HybridNetwork)::Vector{Edge}
    nedges = Edge[];
    for E in N.edge
        E.hybrid && continue
        getchild(E).leaf && continue
        E.length < 0 && push!(nedges, E)
    end
    return nedges
end


"""
    deepcopynetwork(net::HybridNetwork)

Creates a "deep" copy of the network `net`, only
copying objects and values that are relevant to
the SNaQ algorithm as it is implemented here.

WARNING: Saves SIGNIFICANT time and memory over
`Base.deepcopy`, but only duplicates portions that
are relevant to the SNaQ algorithm, so does not
create a true deepcopy.
"""
function deepcopynetwork(net::HybridNetwork)::HybridNetwork
    # List of nodes W/O attached edges
    node_map::Dict{Int, Node} = Dict{Int, Node}()
    sizehint!(node_map, net.numnodes)   # else it rehashes its way up on every copy
    nodec = Array{Node}(undef, net.numnodes)
    for (j, node) in enumerate(net.node)
        nodec[j] = Node(node.number, node.leaf, node.hybrid)
        nodec[j].name = node.name
        node_map[node.number] = nodec[j]
    end

    # List of edges - also attached the edges to their respective nodes
    edgec = Array{Edge}(undef, net.numedges)
    for (j, e) in enumerate(net.edge)
        n1, n2 = node_map[e.node[1].number], node_map[e.node[2].number]
        edgec[j] = Edge(e.number, e.length, e.hybrid, e.gamma, [n1, n2])
        edgec[j].ischild1 = e.ischild1
        edgec[j].ismajor = e.ismajor
        edgec[j].containroot = e.containroot
        push!(n1.edge, edgec[j])
        push!(n2.edge, edgec[j])
    end

    # List of hybrids
    hybc = Array{Node}(undef, net.numhybrids)
    for (j, hyb) in enumerate(net.hybrid)
        hybc[j] = node_map[hyb.number]
    end

    # List of leaves
    leafc = Array{Node}(undef, length(net.leaf))
    for (j, l) in enumerate(net.leaf)
        leafc[j] = node_map[l.number]
    end

    netc = HybridNetwork()
    netc.numtaxa = net.numtaxa
    netc.numnodes = net.numnodes
    netc.numedges = net.numedges
    netc.node = nodec
    netc.edge = edgec
    netc.leaf = leafc
    netc.rooti = net.rooti
    netc.names = net.names
    netc.hybrid = hybc
    netc.numhybrids = net.numhybrids
    netc.isrooted = net.isrooted
    SNaQscore!(netc, SNaQscore(net))
    return netc
end
