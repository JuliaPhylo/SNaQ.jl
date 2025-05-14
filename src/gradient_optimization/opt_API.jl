
const IdxObjMap = Dict{Int, Union{Node, Edge}};   # for readability


function optimize_topology!(
    Nprime::HybridNetwork,
    old_eqns::Vector{QuartetData},
    Nprime_eqns::Vector{QuartetData},
    move::Symbol,
    params::Tuple,
    q::Matrix{Float64},
    q_idxs::Vector{Int},
    opt_maxeval::Int,
    force_resample_all::Bool,
    rng::TaskLocalRNG,
    α::Float64
)

    Nprime_eqns::Vector{QuartetData} = Array{QuartetData}(undef, 3*length(old_eqns))
    if !force_resample_all && can_update_inplace(move)
        @debug "\tGathering updated quartet equations."
        _, param_map, idxobjmap, _ = gather_optimization_info(Nprime, true)
        update_quartet_equations!(old_eqns, Nprime_eqns, Nprime, param_map, move, params, α)
    else
        @debug "\tGathering quartet equations."
        find_quartet_equations!(Nprime, q_idxs, Nprime_eqns);
    end

    @debug "\tOptimizing branch lengths."
    Nprime_eqns = Nprime_eqns[1:length(old_eqns)]
    Nprime_logPL = optimize_bls!(Nprime, Nprime_eqns, q, α; maxeval=opt_maxeval)

    return Nprime_logPL, Nprime_eqns
end


function compute_loss(N::HybridNetwork, q, q_idxs::Vector{Int}, rng::TaskLocalRNG, α::Real)::Tuple{Float64, Vector{QuartetData}}
    @debug "\tGathering quartet equations."
    N_qdata, _, N_params, _ = find_quartet_equations(N, q_idxs)

    @debug "\tComputing loss."
    return compute_loss(N_qdata, N_params, q[q_idxs, :], α), N_qdata
end




function optimize_bls!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::AbstractVector{<:PhyloNetworks.QuartetT},
    α::Real=Inf;
    maxeval::Int=25
)
    obsCF_static = Array{Float64}(undef, length(observed_CFs), 3)
    for j = 1:length(observed_CFs)
        for k = 1:3
            obsCF_static[j, k] = observed_CFs[j].data[k]
        end
    end
    return optimize_bls!(net, eqns, obsCF_static, α, maxeval=maxeval)
end


function optimize_bls!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::Matrix{Float64},
    α::Real=Inf;
    maxeval::Int=10
)

    narg, param_map, idx_obj_map, params, LB, UB, init_steps = gather_optimization_info(net, false)
    #opt = Opt(NLopt.LD_TNEWTON_PRECOND, narg)  # more accurate, but takes longer
    opt = Opt(NLopt.LD_LBFGS, narg)     # faster, but less accurate

    opt.maxeval = maxeval
    opt.ftol_rel = 1e-12
    opt.ftol_abs = 1e-12
    opt.xtol_rel = 1e-8
    opt.xtol_abs = 1e-8

    initial_step!(opt, init_steps)
    opt.lower_bounds = LB
    opt.upper_bounds = UB

    setX!(net, [min(ub / 2.0, val) for (ub, val) in zip(UB, params)], idx_obj_map)
    NLopt.max_objective!(opt, (x, grad) -> objective(x, grad, net, eqns, observed_CFs, idx_obj_map, α))
    (minf, minx, ret) = NLopt.optimize(opt, params)
    setX!(net, minx, idx_obj_map)
    if minf == -Inf
        @show writenewick(net, round=true)
        @info "\n\n\n"
        @show minf
        @show minx
        @show ret
        @show params
        @show UB
        @show LB
        @show UB .- LB
        @show objective(minx, fill(0.0, length(minx)), net, eqns, observed_CFs, idx_obj_map, α)
        error("minf == -Inf")
    end

    # The major/minor property of some hybrid edges may need to be changed at this point
    for hyb in net.hybrid
        par = getparents(hyb)
        if getconnectingedge(hyb, par[1]).gamma > 0.5
            getconnectingedge(hyb, par[1]).ismajor = true
            getconnectingedge(hyb, par[2]).ismajor = false
        elseif getconnectingedge(hyb, par[1]).gamma < 0.5
            getconnectingedge(hyb, par[1]).ismajor = false
            getconnectingedge(hyb, par[2]).ismajor = true
        end

        # If they are both exactly 0.5, just leave the values as they were before.
        # This way, no updates will be forced.
    end

    return minf
end
optimize_bls!(net::HybridNetwork, oCFs; kwargs...) = optimize_bls!(net, find_quartet_equations(net)[1], oCFs; kwargs...)


function objective(X::Vector{T}, grad::Vector{T}, net::HybridNetwork, eqns::Array{QuartetData}, obsCFs::Matrix{T}, idx_obj_map::IdxObjMap, α::Float64)::Float64 where T<:Float64
    setX!(net, X, idx_obj_map)
    fill!(grad, 0.0)
    loss = compute_loss_and_gradient!(eqns, X, grad, obsCFs, α)
    return loss
end



"""
Sets the branch lengths and γ values of edges in `net` according
to the values provided in `X` and the index-to-object map
provided in `idx_obj_map`.
"""
function setX!(net::HybridNetwork, X::Vector{Float64}, idx_obj_map::IdxObjMap)
    for j = 1:length(X)
        obj::Union{Node, Edge} = idx_obj_map[j]
        if typeof(obj) <: PN.Node
            E_major = getparentedge(obj)
            E_minor = getparentedgeminor(obj)
            E_major.gamma = 1-X[j]
            E_minor.gamma = X[j]
        else
            obj.length = X[j]
        end
    end
end



function gather_optimization_info(net::HybridNetwork, change_numbers::Bool=true)

    param_map = Dict{Int, Int}()
    idx_obj_map::IdxObjMap = IdxObjMap();
    uq_ID = net.numedges
    param_idx = 1

    if change_numbers
        for obj in vcat(net.hybrid, net.edge, net.node)
            obj.number = uq_ID
            uq_ID += 1
        end
    end

    order = sortperm([obj.number for obj in vcat(net.hybrid, net.edge)])
    for obj in vcat(net.hybrid, net.edge)[order]
        if typeof(obj) <: Edge && getchild(obj).leaf continue end

        haskey(param_map, obj.number) && error("Duplicate object number #$(obj.number).")
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
            LB[j] = 0.0
            UB[j] = 1.0
            init_steps[j] = 0.1
        else
            params[j] = obj.length
            LB[j] = 0.0
            UB[j] = 25.0
            init_steps[j] = 1.0
        end
    end

    return narg, param_map, idx_obj_map, params, LB, UB, init_steps
end


function gather_params(net::HybridNetwork, param_map::Dict{Int, Int})::Array{Float64}
    params = zeros(length(param_map))
    for obj in vcat(net.hybrid, net.edge)
        if haskey(param_map, obj.number)
            params[param_map[obj.number]] = typeof(obj) <: Node ? getparentedgeminor(obj).gamma : obj.length
        end
    end
    return params
end


function gather_params(net::HybridNetwork)::Array{Float64}
    param_map = gather_optimization_info(net, true)[2]
    return gather_params(net, param_map)
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


function compute_eCFs(net::HybridNetwork, α::Real=Inf)
    eqns, _, params, _ = find_quartet_equations(net)
    eCFs = zeros(length(eqns), 3)
    for j = 1:size(eCFs)[1]
        eCFs[j, 1], eCFs[j, 2] = compute_eCF(eqns[j], params, α)
        eCFs[j, 3] = 1 - eCFs[j, 1] - eCFs[j, 2]
    end
    return eCFs
end