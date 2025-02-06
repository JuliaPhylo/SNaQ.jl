# using Graphs, PhyloNetworks
using PhyloNetworks

# snaq!(tre0, df, restriction_set(max_level=3, require_galled_tree=true))

function restriction_set(; max_level::Real=Inf, require_galled_tree::Bool=false, require_galled_network::Bool=false)

    restrictions = Vector{Function}()
    if max_level < Inf
        push!(restrictions, restrict_maximum_level(max_level))
    end

    if require_galled_tree
        push!(restrictions, restrict_galled_tree())
    elseif require_galled_network
        push!(restrictions, restrict_galled_network())
    end
    return (net) -> all(F(net) for F in restrictions)

end

restrict_maximum_level(level::Int) = (net) -> get_network_level(net) <= level
restrict_galled_network() = (net) -> PhyloNetworks.isgalled(net)
restrict_galled_tree() = (net) -> get_network_level(net) <= 1
restrict_rooted_tree_child() = (net) -> PhyloNetworks.istreechild(net)[1]
restrict_weakly_tree_child() = (net) -> PhyloNetworks.istreechild(net)[2]
restrict_strongly_tree_child() = (net) -> PhyloNetworks.istreechild(net)[3]


function meets_constraints(net::HybridNetwork, constraints::Vector{Function})
    for F in constraints
        if !F(net) return false end
    end
    return true
end


function get_network_level(net::HybridNetwork)
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
is_galled_tree(net::HybridNetwork) = get_network_level(net) <= 1