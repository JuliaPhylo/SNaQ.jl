using NLopt
include("misc.jl")
include("CF_struct.jl")
include("CF_math.jl")


function optimize_bls!(net::HybridNetwork, quartet_eqns, observed_CFs; return_logPL::Bool=false)
    narg = nopt_params(net)
    ngamma = net.numhybrids
    optmap, idxobjmap = generate_optimization_map(net)

    objective(fill(0.5, narg), fill(0.5, narg), net, idxobjmap, quartet_eqns, observed_CFs, optmap)

    opt = Opt(NLopt.LD_LBFGS, narg)
    opt.lower_bounds = [(j <= ngamma) ? 0.0001 : 0.0001 for j=1:narg]
    opt.upper_bounds = [(j <= ngamma) ? 0.5 : 15.0 for j = 1:narg]
    NLopt.max_objective!(opt, (x, grad) -> objective(x, grad, net, idxobjmap, quartet_eqns, observed_CFs, optmap))
    (minf, minx, ret) = NLopt.optimize(opt, fill(0.5, narg))

    if return_logPL
        return get_params(net), compute_logPL(quartet_eqns, net.edge, observed_CFs)
    end
    return get_params(net)
end
optimize_bls!(net::HybridNetwork, oCFs; kwargs...) = optimize_bls!(net, find_quartet_equations(net)[1], oCFs; kwargs...)


function objective(X::Vector{Float64}, grad::Vector{Float64}, N, idx_to_obj_map, q_eqns, obsCFs, opt_map)::Float64
    for j = 1:length(X)
        obj::Union{PN.Edge, PN.Node} = idx_to_obj_map[j]
        if typeof(obj) <: PN.Node
            getparentedge(obj).gamma = 1-X[j]
            getparentedgeminor(obj).gamma = X[j]
        else
            obj.length = X[j]
        end
    end

    total_loss = 0.0
    for ((eqns1, eqns2, eqns3), iter_oCF) in zip(q_eqns, obsCFs)
        total_loss += iter_oCF.data[1] * log(compute_eCF(eqns1, N.edge))
        total_loss += iter_oCF.data[2] * log(compute_eCF(eqns2, N.edge))
        total_loss += iter_oCF.data[3] * log(compute_eCF(eqns3, N.edge))
    end
    grad .= compute_gradient(q_eqns, N.edge, opt_map, obsCFs)
    return total_loss
end

