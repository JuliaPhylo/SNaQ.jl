"""
Generates a `Dict` mapping all edges and hybrid nodes w/ an associated γ parameter to their
index in the vector that is being optimized.
"""
function generate_optimization_map(net::HybridNetwork)::Tuple{Dict{Union{PN.Edge,PN.Node}, Int}, Dict{Int, Union{PN.Edge,PN.Node}}}
    ngammas = length(net.hybrid)
    nedges = net.numedges - net.numtaxa
    map = Dict{Union{PN.Edge,PN.Node}, Int}()
    j = 1
    for H in net.hybrid
        map[H] = j
        j += 1
    end
    for E in net.edge
        getchild(E).leaf && continue
        # E.hybrid && !E.ismajor && continue
        map[E] = j
        j += 1
    end
    return map, Dict{Int, Union{PN.Edge,PN.Node}}(map[key] => key for key in keys(map))
end
nopt_params(net::HybridNetwork)::Int = length(generate_optimization_map(net)[1])


"""
Function mostly used for debugging, gets the relevant parameters for network `net`.
"""
function get_params(net::HybridNetwork)
    vals = []
    optmap, _ = generate_optimization_map(net)
    for key in keys(optmap)
        if typeof(key) <: PN.Node
            push!(vals, getparentedgeminor(key).gamma)
        else
            push!(vals, key.length)
        end
    end
    return vals
end