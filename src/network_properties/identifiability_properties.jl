
"""
Checks whether network `N` has an identifiable semi-directed topology
and whether all of its branch parameters are known to be identifiable.
Returns `true` if the topology and branch parameters are BOTH identifiable,
`false` if either or both are not identifiable.

Checklist:
- Rooted --> root edge lengths are non-identifiable
- All trees --> identifiable
- 2-cycles --> non-identifiable topology
- 3-cycles --> non-identifiable topology
- Level-1:
    - Bad diamonds --> non-identifiable branch lengths
- Tree-child galled:
    - Network is C5 --> fully identifiable
    - Network is C4 --> fully identifiable if and only if 2 samples per taxa
    - Otherwise --> unknown, return false
- Non tree-child galled --> unknown, return false
"""
function knownidentifiable(N::HybridNetwork)::Bool
    !N.isrooted || return false
    N.numhybrids == 0 && return true

    N = deepcopy_network(N)
    (shrink2cycles!(N) || shrink3cycles!(N)) && return false

    if getlevel(N) == 1
        # Only need to check for bad diamonds in each blob
        # Diamonds are GOOD if and only if
        # (n0 >= 2 OR n2 >= 2) OR (n1 >= 2 AND n3 >= 2)
        # (see https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896#pgen-1005896-g005)
        blobs = [bcp for bcp in biconnectedcomponents(N) if length(bcp) > 1];
        for blob in blobs
            cycle_nodes = unique(reduce(vcat, [[getchild(e), getparent(e)] for e in blob]))
            length(cycle_nodes) == 4 || continue    # if this is not a diamond, it cannot be a 3-cycle,
                                                    # so must have at least 5 nodes, and so is fine
            n0, n1, n2, n3 = getdiamondns(cycle_nodes)
            if !(n0 >= 2 || n2 >= 2 || (n1 >= 2 && n3 >= 2))
                return false
            end
        end
    else
        boolCk::Bool, k::Int = isCk(N)
        boolCk || return false
        k < 4 && return false
        k > 4 && return true

        # if network is C4 we need 2 samples per taxon,
        # ie. at least 2 taxa under each lowest nodes
        # in each blob. We know that the network is C4
        # at this point, so we don't have to bother checking
        # every blob, all of the lowest nodes are exactly
        # the set of the hybrids in `N`.
        lnodes::Vector{Node} = N.hybrid
        for ln in lnodes
            if length(getleavesbelow([ln])) < 2
                return false
            end
        end
    end
    return true
end


"""
Gets the `n0,n1,n2,n3` taxa counts associated with the 4-cycle `cycle`.
Order of `cycle` does not matter. Counts returned in the order `n0,n1,n2,n3`.
"""
function getdiamondns(cycle::Vector{Node})::Vector{Int}
    length(cycle) == 4 || error("cycle must be of length 4 (has length $(length(cycle)))")
    
    n3root::Node = cycle[findfirst(node -> node.hybrid, cycle)]
    n1root::Node = cycle[findfirst(node -> node != n3root && all(e -> e.node[1] != n3root && e.node[2] != n3root, node.edge), cycle)]
    n0root::Node, n2root::Node = cycle[findall(node -> node != n3root && node != n1root && node != n1root && node != n3root, cycle)]

    ns::Vector{Int} = [0, 0, 0, 0]
    for (nidx, rootnode) in enumerate([n0root, n1root, n2root, n3root])
        visitednodes::Vector{Node} = [n0root, n1root, n2root, n3root]
        Q::Vector{Node} = [rootnode]
        
        while length(Q) > 0
            currnode::Node = Q[1]
            deleteat!(Q, 1)
            
            push!(visitednodes, currnode)
            if currnode.leaf
                ns[nidx] += 1
            else
                for anode in attatchednodes(currnode)
                    if anode ∈ visitednodes continue end
                    push!(Q, anode)
                end
            end
        end
    end
    return ns
end


attatchednodes(n::Node) = [e.node[1] == n ? e.node[2] : e.node[1] for e in n.edge]


"""
Returns `(false, 0)` if network `N` is not of any class Ck and
    returns `(true, k)` otherwise.
"""
function isCk(N::HybridNetwork)::Tuple{Bool, Int}
    N = deepcopy_network(N)
    is_galled_network(N) || return (false, 0)
    istreechild(N)[3] || return (false, 0)
    if shrink2cycles!(N) return (false, 0) end

    k::Int = N.numtaxa
    blobs = [bcp for bcp in biconnectedcomponents(N) if length(bcp) > 1]
    for blob in blobs
        # lnodes are the "lowest nodes in the blob B". Because we know at this
        # point that N is galled, these nodes must be exactly the set of
        # hybrids in the blob.
        lnodes = unique([getchild(e) for e in blob if getchild(e).hybrid])
        Y_l = getleavesbelow(lnodes)
        for l = 1:length(Y_l)
            bloblet = deepcopy_network(N)
            for i = 1:length(Y_l)
                if i == l continue end
                deleteleaf!(bloblet, Y_l[i].name; simplify=true, nofuse=false, multgammas=false)
            end
            bloblet_bcp = [bcp for bcp in biconnectedcomponents(bloblet) if length(bcp) > 1]
            if length(bloblet_bcp) != 1
                return (false, 0)
            end
            k = min(k, length(bloblet_bcp[1]))
        end
    end
    return (true, k)
end


"""
Checks whether network `N` is of class Ck given int `k`.
"""
function isCk(N::HybridNetwork, k::Int)::Bool
    class::Tuple{Bool, Int} = isCk(N)
    return class[1] && class[2] >= k
end


"""
Gets all of the leaves BELOW each of the nodes in the vector `nodes`.
Even if the network is semi-directed, traverses the edges as though
they are directed.
"""
function getleavesbelow(nodes::Vector{Node})::Vector{Node}
    nodes::Vector{Node} = [node for node in nodes]
    leaves::Vector{Node} = []
    nodes_visited::Vector{Node} = []

    while length(nodes) > 0
        curr_node::Node = nodes[1]
        deleteat!(nodes, 1)
        curr_node in nodes_visited && continue
        push!(nodes_visited, curr_node)
        
        if curr_node.leaf
            push!(leaves, curr_node)
            continue
        end

        for n in getchildren(curr_node)
            push!(nodes, n)
        end
    end
    return leaves
end



# # c3 is NOT identifiable
# c3 = readnewick("((c,#H1),((((g)#H1,(d,#H4)),#H3),(#H2,((b2)#H4,((f)#H2,((b1)#H3,a))))));")
# SNaQ.semidirect_network!(c3)
# @test !SNaQ.knownidentifiable(c3)

# # c4 is NOT identifiable
# c4 = deepcopy(not_c4); deleteleaf!(c4, "f")
# SNaQ.semidirect_network!(c4)
# @test !SNaQ.knownidentifiable(c4)

# # c6 is identifiable
# c6 = readnewick("(((((((a1,a2))#H1,g),f),#H2),e),(((((b1,b2))#H2,c),#H1),d));")
# SNaQ.semidirect_network!(c6)
# @test SNaQ.knownidentifiable(c6)

# # dianet is identifiable
# dianet = readnewick("((a,(b,#H1)),(((c,d))#H1,(e,f)));");
# SNaQ.semidirect_network!(dianet);
# @test SNaQ.knownidentifiable(dianet)