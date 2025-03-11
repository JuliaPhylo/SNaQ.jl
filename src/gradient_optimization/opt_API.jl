using NLopt
include("misc.jl")
include("CF_struct.jl")
include("CF_math.jl")




function optimize_bls!(net::HybridNetwork, quartet_eqns, observed_CFs; return_logPL::Bool=false, α::Real=Inf)

    narg = length(net.hybrid) + length(net.edge)
    param_map = Dict{Int, Int}()
    idx_obj_map = Dict{Int, Union{Node, Edge}}()
    params = Array{Float64}(undef, narg)
    LB = Array{Float64}(undef, narg)
    UB = Array{Float64}(undef, narg)
    max_ID = net.numedges

    for (j, obj) in enumerate(vcat(net.hybrid, net.edge))
        obj.number = max_ID + j
        param_map[obj.number] = j
        idx_obj_map[j] = obj

        if typeof(obj) <: Node
            params[j] = getparentedgeminor(obj).gamma
            LB[j] = 0.0001
            UB[j] = 0.9999
        else
            params[j] = obj.length
            LB[j] = 0.0001
            UB[j] = 15.0
        end
    end

    blocks, _ = find_quartet_equations(net)

    opt = Opt(NLopt.LD_LBFGS, narg)
    opt.lower_bounds = LB
    opt.upper_bounds = UB
    NLopt.max_objective!(opt, (x, grad) -> objective(x, grad, net, blocks, observed_CFs, idx_obj_map, α))
    (minf, minx, ret) = NLopt.optimize(opt, fill(0.5, narg))

    return minx
end
optimize_bls!(net::HybridNetwork, oCFs; kwargs...) = optimize_bls!(net, find_quartet_equations(net)[1], oCFs; kwargs...)


function objective(X::Vector{Float64}, grad::Vector{Float64}, net::HybridNetwork, blocks::AbstractArray, obsCFs, idx_obj_map, α)::Float64
    for j = 1:length(X)
        obj::Union{Node, Edge} = idx_obj_map[j]
        if typeof(obj) <: PN.Node
            E_major = getparentedge(obj)
            E_minor = getparentedgeminor(obj)
            E_major.gamma = 1-X[j]
            E_minor.gamma = X[j]

            if 1-X[j] < 0.5
                E_major.ismajor = false
                E_minor.ismajor = true
            end
        else
            obj.length = X[j]
        end
    end

    total_loss = 0.0
    iter_eCF::Float64 = 0.0
    for j = 1:size(blocks)[1]
        for k = 1:3
            # Loss function
            iter_eCF = compute_eCF(blocks[j,k], X, k, α)
            total_loss += obsCFs[j].data[k] * log(iter_eCF)

            # Gradient
            block_derivs = compute_block_derivs(blocks[j,k], X, k, α)
            grad .+= block_derivs .* block_derivs .* obsCFs[j].data[k] ./ iter_eCF
        end
    end
    return total_loss
end

