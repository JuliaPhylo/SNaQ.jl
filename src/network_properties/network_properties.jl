using PhyloNetworks

# snaq!(tre0, df, restrictionset(max_level=3, require_galled_tree=true))

function restrictionset(; max_level::Real=Inf, galled_tree::Bool=false, galled_network::Bool=false,
    rooted_tree_child::Bool=false, weakly_tree_child::Bool=false, strongly_tree_child::Bool=false)

    restrictions = Vector{Function}()
    if max_level < Inf
        push!(restrictions, restrictmaximumlevel(max_level))
    end

    if galled_tree
        push!(restrictions, restrictgalledtree())
    elseif galled_network
        push!(restrictions, restrictgallednetwork())
    end

    if rooted_tree_child
        push!(restrictions, restrictrootedtreechild())
    elseif weakly_tree_child
        push!(restrictions, restrictweaklytreechild())
    elseif strongly_tree_child
        push!(restrictions, restrictstronglytreechild())
    end

    return (net) -> all(F(net) for F in restrictions)

end

restrictmaximumlevel(level::Int) = (net) -> getnetworklevel(net) <= level
restrictgallednetwork() = (net) -> isgallednetwork(net)
restrictgalledtree() = (net) -> getnetworklevel(net) <= 1
restrictrootedtreechild() = (net) -> PhyloNetworks.istreechild(net)[1]
restrictweaklytreechild() = (net) -> PhyloNetworks.istreechild(net)[2]
restrictstronglytreechild() = (net) -> PhyloNetworks.istreechild(net)[3]
defaultrestrictions() = (net) -> knownidentifiable(net)
norestrictions() = (net) -> true


function meets_constraints(net::HybridNetwork, constraints::Vector{Function})
    for F in constraints
        if !F(net) return false end
    end
    return true
end


# Faster than PhyloNetworks but not as extensively tested, so we don't use it
function getnetworklevel(net::HybridNetwork)
    bi_comp = [comp for comp in biconnectedcomponents(net) if length(comp) > 1]
    max_level = 1
    for component in bi_comp
        hybrids_seen = [];
        for edge in component
            child = getchild(edge)
            if child.hybrid && !(child in hybrids_seen)
                push!(hybrids_seen, child)
            end
        end
        max_level = max(max_level, length(hybrids_seen))
    end
    return max_level
end
isgalledtree(net::HybridNetwork) = getnetworklevel(net) <= 1


function isgallednetwork(net::HybridNetwork)
    bi_comp = [comp for comp in biconnectedcomponents(net) if length(comp) > 1]
    v_counts = Dict{Node, Int}()
    for component in bi_comp
        component = reduce(vcat, [[edge.node[1], edge.node[2]] for edge in component])
        for node in component
            if node.hybrid
                if !haskey(v_counts, node)
                    v_counts[node] = 1
                else
                    v_counts[node] += 1
                    if v_counts[node] > 2 return false end
                end
            end
        end
    end
    return true
end
