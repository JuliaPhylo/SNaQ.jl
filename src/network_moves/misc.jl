using PhyloNetworks;


function semidirect_network!(N::HybridNetwork; check_conditions::Bool=true)
    # Semi-directed network attribute that we are enforcing here:
    #   all non-leaf vertices have degree-3
    length(getroot(N).edge) == 2 || error("Expected the root to have 2 attached edges.")
    sum([E.hybrid for E in getroot(N).edge]) <= 1 || error("Both root edges are hybrids.")

    # Fuse root edges and mark the network as no longer rooted
    PhyloNetworks.fuseedgesat!(N.rooti, N, false)
    N.isrooted = false

    if check_conditions
        # Check that all network assumptions are met
        for n in N.node
            if n.leaf
                length(n.edge) == 1 || error("Found leaf in `N` with degree $(length(n.edge))")
            else
                length(n.edge) == 3 || error("Found internal node in `N` with degrees $(length(n.edge))")
            end
        end
        all(length(e.node) == 2 for e in N.edge) || error("Found edge with ≠2 attached nodes.")
    end
end


"""
Helper function that swaps all hybridization-related information between
    edges `e1` and `e2`. Used heavily in the network move functions to
    ensure that manipulated edges have valid reticulation info.
    `e1` and `e2` needn't be reticulations.
"""
function swap_edge_hybrid_info!(e1::Edge, e2::Edge)
    temp_e2_gamma = e2.gamma
    temp_e2_hybrid = e2.hybrid
    temp_e2_ismajor = e2.ismajor

    e2.gamma = e1.gamma
    e2.hybrid = e1.hybrid
    e2.ismajor = e1.ismajor

    e1.gamma = temp_e2_gamma
    e1.hybrid = temp_e2_hybrid
    e1.ismajor = temp_e2_ismajor
end