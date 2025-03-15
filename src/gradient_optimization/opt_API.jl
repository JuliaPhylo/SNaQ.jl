using NLopt
include("misc.jl")
include("CF_struct.jl")
include("CF_blocks.jl")
include("CF_recursive_blocks.jl")
include("CF_equations.jl")




function optimize_bls!(net::HybridNetwork, blocks, observed_CFs, α::Real=Inf; return_logPL::Bool=false)

    narg, param_map, idx_obj_map, params, LB, UB, init_steps = gather_optimization_info(net)

    opt = Opt(NLopt.LD_LBFGS, narg)

    opt.maxeval = 25
    opt.ftol_rel = 1e-6
    opt.ftol_abs = 1e-6
    opt.xtol_rel = 1e-3
    opt.xtol_abs = 1e-3

    opt.initial_step = init_steps
    opt.lower_bounds = LB
    opt.upper_bounds = UB

    NLopt.max_objective!(opt, (x, grad) -> objective(x, grad, net, blocks, observed_CFs, idx_obj_map, α))
    (minf, minx, ret) = NLopt.optimize(opt, fill(0.1, narg))

    return minf
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

    fill!(grad, 0.0)
    total_loss::Float64 = compute_loss_and_gradient!(blocks, X, grad, obsCFs, α)


    # total_loss::Float64 = 0.0
    # iter_eCF::Float64 = 0.0
    # iter_derivs::Array{Float64} = zeros(length(X))
    # for j = 1:size(blocks)[1]
    #     for k = 1:3
    #         fill!(iter_derivs, 0.0)

    #         # Loss function
    #         iter_eCF = compute_block_eCF_and_gradient!(blocks[j, k], X, iter_derivs, k, α)
    #         iter_eCF = max(iter_eCF, 1e-9)

    #         if obsCFs[j].data[k] > 0
    #             # if obsCF is about 0, skip this (i.e. add 0)
    #             total_loss += obsCFs[j].data[k] * log(iter_eCF / obsCFs[j].data[k])
    #         end

    #         # Gradient
    #         grad .+= obsCFs[j].data[k] .* iter_derivs ./ iter_eCF
    #     end
    # end
    return total_loss
end


function gather_optimization_info(net::HybridNetwork)

    param_map = Dict{Int, Int}()
    idx_obj_map = Dict{Int, Union{Node, Edge}}()
    max_ID = net.numedges
    param_idx = 1

    for obj in vcat(net.hybrid, net.edge)
        if typeof(obj) <: Edge && getchild(obj).leaf continue end

        obj.number = max_ID + param_idx
        param_map[obj.number] = param_idx
        idx_obj_map[param_idx] = obj
        param_idx += 1
    end
    
    params = gather_params(net, param_map)
    narg = length(param_map)
    LB = Array{Float64}(undef, narg)
    UB = Array{Float64}(undef, narg)
    init_steps = Array{Float64}(undef, narg)

    for j = 1:narg
        obj::Union{Node,Edge} = idx_obj_map[j]
        if typeof(obj) <: Node
            params[j] = getparentedgeminor(obj).gamma
            LB[j] = 0.0001
            UB[j] = 0.9999
            init_steps[j] = 0.1
        else
            params[j] = obj.length
            LB[j] = 0.0
            UB[j] = 15.0
            init_steps[j] = 1.0
        end
    end

    return narg, param_map, idx_obj_map, params, LB, UB, init_steps
end


function gather_params(net::HybridNetwork, param_map::Dict{Int, Int})::Array{Float64}
    params = zeros(length(param_map))
    param_idx = 1
    for obj in vcat(net.hybrid, net.edge)
        if haskey(param_map, obj.number)
            params[param_idx] = typeof(obj) <: Node ? getparentedgeminor(obj).gamma : obj.length
        end
    end
    return params
end


function compute_gradient(net::HybridNetwork, obsCFs)
    blocks, _, params, _ = find_quartet_equations(net)
    
    grad = zeros(length(params))
    for j = 1:size(blocks)[1]
        for k = 1:3
            iter_eCF = compute_eCF(blocks[j, k], params, k, Inf)
            block_derivs = compute_block_derivs(blocks[j, k], params, k, Inf)
            grad .+= obsCFs[j].data[k] .* block_derivs ./ iter_eCF
        end
    end
    return grad
end


function compute_logPL(blocks::AbstractMatrix{<:AbstractArray{<:AbstractArray{<:Block}}}, params::AbstractArray{<:Real}, obsCFs, α::Real)::Float64
    total_loss::Float64 = 0.0
    for j = 1:size(blocks)[1]
        for k = 1:3
            # Loss function
            total_loss += obsCFs[j].data[k] * log(compute_eCF(blocks[j,k], params, k, α) / obsCFs[j].data[k])
        end
    end
    return total_loss
end


function compute_eCFs(net::HybridNetwork)
    blocks, _, params, _ = find_quartet_equations(net)
    eCFs = zeros(size(blocks))
    for j = 1:size(eCFs)[1]
        for k = 1:3
            eCFs[j, k] = compute_eCF(blocks[j, k], params, k, Inf)
        end
    end
    return eCFs
end