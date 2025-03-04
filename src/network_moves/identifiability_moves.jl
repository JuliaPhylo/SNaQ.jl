using PhyloNetworks
import PhyloNetworks: getOtherNode


function shrink_bad_diamonds!(N::HybridNetwork)

    check_again::Bool = true
    while check_again
        check_again = false
        for H in N.hybrid
            # 1. check if H is apart of a 4-cycle
            p_major = getparent(getparentedge(H))
            p_minor = getparent(getparentedgeminor(H))
            pmaj1, pmaj2 = [getOtherNode(edge, p_major) for edge in p_major.edge if getOtherNode(edge, p_major) != H]
            pmin1, pmin2 = [getOtherNode(edge, p_minor) for edge in p_minor.edge if getOtherNode(edge, p_minor) != H]
            fourth_node = nothing

            if pmin1 == pmaj1 || pmin1 == pmaj2
                fourth_node = pmin1
            elseif pmin2 == pmaj1 || pmin2 == pmaj2
                fourth_node = pmin2
            end

            fourth_node !== nothing || continue

            # 2. check if n1>=2 and n3>=2
            n3_ge2 = length(getchildren(H)) > 1 || !getchild(H).leaf
            if n3_ge2
                n1_subtree_node = get_other_nodes(fourth_node, p_major, p_minor)[1]
                n1_ge2 = !n1_subtree_node.leaf && length(get_other_nodes(n1_subtree_node, fourth_node)) > 1
                n1_ge2 && continue  # this is a good diamond
            end

            # 3. check (2) failed, now check if n1>=2 OR n0>=2
            n2_subtree_node = get_other_nodes(p_major, fourth_node, H)[1]
            n2_ge2 = !n2_subtree_node.leaf && length(get_other_nodes(n2_subtree_node, p_major)) > 1
            n2_ge2 && continue  # good diamond
            
            n0_subtree_node = get_other_nodes(p_minor, fourth_node, H)[1]
            n0_ge2 = !n0_subtree_node.leaf && length(get_other_nodes(n0_subtree_node, p_minor)) > 1
            n0_ge2 && continue  # good diamond

            # 4. this is a bad diamond - shrink it! (i.e. remove the minor reticulate edge)
            remove_hybrid!(H, N)
            check_again = true      # this could drastically change the network, so check everything again!
            break
        end
    end

end


function get_other_nodes(n::Node, except_nodes::Node...)
    return [getOtherNode(e, n) for e in n.edge if !(getOtherNode(e, n) in except_nodes)]
end