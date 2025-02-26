using PhyloNetworks
import PhyloNetworks: breakedge!, fuseedgesat!
include("misc.jl")
include("rNNI_validity.jl")


"""
Performs an rSPR move according to the procedure defined [here in Figure 6](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005611).
"""
function perform_rSPR!(N::HybridNetwork, w::Node, x::Node, y::Node, z::Node, xprime::Node, yprime::Node)
    is_valid_rSPR(w, x, y, z, xprime, yprime) || error("Topological conditions for rSPR not met.")

    which_case = (z.hybrid) ? 2 : 1 # 1 corresponds to the upper half of Figure 6
                                    # 2 corresponds to the lower half of Figure 6

    # 1. Split the x'y' arc, adding a new node z' (just z in Figure 6)
    edge_xprime_yprime = xprime.edge[findfirst(e -> yprime in e.node, xprime.edge)]
    zprime, edge_xprime_zprime = breakedge!(edge_xprime_yprime, N)
    zprime.name = z.name    # maintain original node naming - will be useful for debugging

    # 2. Disconnect z from w and attach w to z'
    edge_zw = z.edge[findfirst(z_edge -> w in z_edge.node, z.edge)]
    replace!(edge_zw.node, z => zprime)
    deleteat!(z.edge, findfirst(z_edge -> z_edge == edge_zw, z.edge))
    push!(zprime.edge, edge_zw)

    # 3. If this is case 2, we need to swap hybrid info before fusing
    if which_case == 2
        edge_xz = z.edge[findfirst(z_edge -> x in z_edge.node, z.edge)]
        swap_edge_hybrid_info!(edge_xprime_zprime, edge_xz)
        zprime.hybrid = true
    end

    # 4. Fuse across the now redundant node z
    z.hybrid = false
    z_idx = findfirst(j -> N.node[j] == z, 1:length(N.node))
    @info N.node[z_idx]
    fuseedgesat!(z_idx, N)
    @info z.edge
end


