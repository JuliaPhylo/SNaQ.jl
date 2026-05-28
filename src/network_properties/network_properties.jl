using PhyloNetworks

# snaq!(tre0, df, restrictionset(max_level=3, require_galled_tree=true))

"""
    restrictionset(; max_level=Inf, galled_tree=false, galled_network=false,
                     rooted_tree_child=false, weakly_tree_child=false, strongly_tree_child=false)

Combine multiple built-in restrictions into a single restriction function for use with
[`snaq!`](@ref). All specified restrictions must be satisfied simultaneously.

If both `galled_tree` and `galled_network` are `true`, `galled_tree` takes precedence.
Similarly, if multiple tree-child options are `true`, only the first one that applies is used:
`rooted_tree_child` is checked first, then `weakly_tree_child`, then `strongly_tree_child`.

# Optional Named Arguments
- `max_level` (default `Inf`): maximum network level (number of hybrids in the most complex biconnected component)
- `galled_tree` (default `false`): restrict to level-1 networks (galled trees)
- `galled_network` (default `false`): restrict to galled networks
- `rooted_tree_child` (default `false`): restrict to rooted tree-child networks
- `weakly_tree_child` (default `false`): restrict to weakly tree-child networks
- `strongly_tree_child` (default `false`): restrict to strongly tree-child networks
"""
function restrictionset(; max_level::Real=Inf, galled_tree::Bool=false, galled_network::Bool=false,
    rooted_tree_child::Bool=false, weakly_tree_child::Bool=false, strongly_tree_child::Bool=false)

    restrictions = Vector{Function}()
    if max_level < Inf
        push!(restrictions, restrictmaximumlevel(max_level))
    end

    if galled_tree
        push!(restrictions, restrictgalledtree)
    elseif galled_network
        push!(restrictions, restrictgallednetwork)
    end

    if rooted_tree_child
        push!(restrictions, restrictrootedtreechild)
    elseif weakly_tree_child
        push!(restrictions, restrictweaklytreechild)
    elseif strongly_tree_child
        push!(restrictions, restrictstronglytreechild)
    end

    return (net) -> all(F(net) for F in restrictions)

end

"""
    restrictmaximumlevel(level)

Return a restriction function that accepts only networks whose level is at most `level`.
The level of a network is the maximum number of hybrid nodes in any biconnected component.
"""
restrictmaximumlevel(level::Int) = (net) -> getnetworklevel(net) <= level

"""
    restrictgallednetwork(net)

Return `true` if `net` is a galled network: each hybrid node appears in at most 2
biconnected components. This is less restrictive than [`restrictgalledtree`](@ref).
"""
restrictgallednetwork(net::HybridNetwork) = isgallednetwork(net)
restrictgallednetwork() = (net) -> isgallednetwork(net)

"""
    restrictgalledtree(net)

Return `true` if `net` is a level-1 network (galled tree): no biconnected component
contains more than one hybrid node.
"""
restrictgalledtree(net::HybridNetwork) = getnetworklevel(net) <= 1
restrictgalledtree() = (net) -> getnetworklevel(net) <= 1

"""
    restrictrootedtreechild(net)

Return `true` if `net` is a rooted tree-child network: every non-leaf node has at least
one child that is not a hybrid node.
"""
restrictrootedtreechild(net::HybridNetwork) = PhyloNetworks.istreechild(net)[1]
restrictrootedtreechild() = (net) -> PhyloNetworks.istreechild(net)[1]

"""
    restrictweaklytreechild(net)

Return `true` if `net` satisfies the weakly tree-child property.
This is less restrictive than [`restrictrootedtreechild`](@ref).
"""
restrictweaklytreechild(net::HybridNetwork) = PhyloNetworks.istreechild(net)[2]
restrictweaklytreechild() = (net) -> PhyloNetworks.istreechild(net)[2]

"""
    restrictstronglytreechild(net)

Return `true` if `net` satisfies the strongly tree-child property.
This is more restrictive than [`restrictrootedtreechild`](@ref).
"""
restrictstronglytreechild(net::HybridNetwork) = PhyloNetworks.istreechild(net)[3]
restrictstronglytreechild() = (net) -> PhyloNetworks.istreechild(net)[3]

"""
    defaultrestrictions(net)

No restrictions applied — equivalent to [`norestrictions`](@ref).
"""
defaultrestrictions(net::HybridNetwork) = norestrictions(net)
defaultrestrictions() = (net) -> norestrictions(net)

"""
    norestrictions(net)

No restrictions applied. Accepts any binary, semi-directed network.
This is the default behavior of [`snaq!`](@ref).
"""
norestrictions(net::HybridNetwork) = true
norestrictions() = (net) -> true


function meetsconstraints(net::HybridNetwork, constraints::Vector{Function})
    for F in constraints
        if !F(net) return false end
    end
    return true
end


# Faster than PhyloNetworks but not as extensively tested, so we don't use it
function getnetworklevel(net::HybridNetwork)
    if net.numhybrids == 0 return 0 end
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
