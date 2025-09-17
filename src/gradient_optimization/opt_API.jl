
const IdxObjMap = Dict{Int, Union{Node, Edge}};   # for readability


"""
Optimizes the branch lengths and γ parameters of the network `Nprime`.

# Arguments:
- `Nprime::HybridNetwork`: network to be optimized, typically the proposed network in [`search`](@ref).
- `old_eqns::Vector{QuartetData}`: quartet equations of the network that came prior to `Nprime` in the
    hill climbing optimization of [`search`](@ref).
- `move::Symbol`: the topological move that turned the previous network into `Nprime`. Used to inform
    which quartets equations in `old_eqns` need to be re-calculated and which remain the same.
- `params::Tuple`: edges and/or hybrid nodes that define the move corresponding to `move` - used along with
    `move` for the aforementioned purpose.
- `q::Matrix{Float64}`: the set of quartet concordance factors used to optimize the network's parameters.
    Should be a matrix with the same number of rows as the length of `old_eqns` and 3 columns.
- `q_idxs::Vector{Int}`: indices corresponding to the set of quartets that are being used for optimization.
    I.e., when `propQuartets` in [`search`](@ref) is 1.0, this contains all integers from 1 to [`nchoose4taxa_length`](@ref).
    When `propQuartets` is 0.1, it contains 10% as many integers, randomly selected in this range.
- `opt_maxeval::Int`: maximum number of evaluations when optimizing parameters.
- `force_resample_all::Bool`: if `true`, ignore `move` and `params` and re-calculate *every* quartet CF
    equation. Typically set to `true` when, e.g., a new reticulation is added to the network.
- `rng::TaskLocalRNG`: `TaskLocalRNG` object from which random numbers are generated. Ensures reproducibility.
- `α::Float64`: inheritance correlation parameter used in calculating pseudo-likelihoods.
"""
function optimize_topology!(
    Nprime::HybridNetwork,
    old_eqns::Vector{QuartetData},
    move::Symbol,
    params::Tuple,
    q::Matrix{Float64},
    q_idxs::Vector{Int},
    opt_maxeval::Int,
    force_resample_all::Bool,
    rng::TaskLocalRNG,
    α::Float64;
    optargs...
)::Tuple{Float64, QuartetData}
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
    Nprime_logPL = optimize_bls!(Nprime, Nprime_eqns, q, α; maxeval=opt_maxeval, optargs...)

    return Nprime_logPL, Nprime_eqns
end


"""
Computes the loss of network `N` given quartet equations `q` sampled from indices `q_idxs`.

Returns a tuple `(loss::Float64, qeqns::Vector{QuartetData})` where `loss` is the aforementioned loss
of network `N` and `qeqns` are the CF equations the define network `N`.
"""
function compute_loss(N::HybridNetwork, q::Matrix{Float64}, q_idxs::Vector{Int}, rng::TaskLocalRNG, α::Real)::Tuple{Float64, Vector{QuartetData}}
    @debug "\tGathering quartet equations."
    N_qdata, _, N_params, _ = find_quartet_equations(N, q_idxs)

    @debug "\tComputing loss."
    return compute_loss(N_qdata, N_params, q[q_idxs, :], α), N_qdata
end



"""
Optimizes the branch lengths of network `net` which is defined by quartet concordance factor
equations `eqns` based on observed quartet CFs `observed_CFs` under inheritance correlation
parameter `α`. `maxeval` adjusts the maximum number of loss evaluations in the optimization
process. This overloaded function is used as a helper method that can directly take 
Vectors of `PhyloNetworks.QuartetT` as input so that the user doesn't need to handle
converting the input data.
"""
function optimize_bls!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::AbstractVector{<:PhyloNetworks.QuartetT},
    α::Real=Inf;
    maxeval::Int=25
)::Float64
    obsCF_static = Array{Float64}(undef, length(observed_CFs), 3)
    for j in eachindex(observed_CFs)
        for k = 1:3
            obsCF_static[j, k] = observed_CFs[j].data[k]
        end
    end
    return optimize_bls!(net, eqns, obsCF_static, α, maxeval=maxeval)
end


"""
Optimizes the branch lengths (and γ parameters) of network `net`.

# Required Arguments:
- `net::HybridNetwork`: the network to be optimized
- `eqns::Array{QuartetData}`: the set of quartet CF equations the define `net`
- `observed_CFs::Matrix{Float64}`: the corresponding set of observed CFs. Each row
    in this matrix MUST line up with each entry in `eqns`. If this is called through
    the [`search`](@ref) method, then this is handled.

# Optional Arguments:
- `α::Float64=Inf`: the inheritance correlation parameter with which the network is optimized. (Default=`Inf`)
- `maxeval::Int=10`: the maximum number of loss evaluations during optimization. (Default=`10`)
- `ftolRel::Float64=1e-12`: optimization parameter passed to the NLOpt.jl optimizer. (Default=`1e-12`)
- `ftolAbs::Float64=1e-12`: optimization parameter passed to the NLOpt.jl optimizer. (Default=`1e-12`)
- `xtolRel::Float64=1e-8`: optimization parameter passed to the NLOpt.jl optimizer. (Default=`1e-8`)
- `xtolAbs::Float64=1e-8`: optimization parameter passed to the NLOpt.jl optimizer. (Default=`1e-8`)
"""
function optimize_bls!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::Matrix{Float64},
    α::Real=Inf;
    maxeval::Int=10,
    ftolRel::Float64=1e-12,
    ftolAbs::Float64=1e-12,
    xtolRel::Float64=1e-8,
    xtolAbs::Float64=1e-8
)::Float64

    narg, param_map, idx_obj_map, params, LB, UB, init_steps = gather_optimization_info(net, false)
    #opt = Opt(NLopt.LD_TNEWTON_PRECOND, narg)  # more accurate, but takes longer
    opt = Opt(NLopt.LD_LBFGS, narg)     # faster, but less accurate

    opt.maxeval = maxeval
    opt.ftol_rel = ftolRel
    opt.ftol_abs = ftolAbs
    opt.xtol_rel = xtolRel
    opt.xtol_abs = xtolAbs

    initial_step!(opt, init_steps)
    opt.lower_bounds = LB
    opt.upper_bounds = UB

    x0::Vector{Float64} = [min(ub / 2.0, val) for (ub, val) in zip(UB, params)]
    NLopt.max_objective!(opt, (x, grad) -> objective(x, grad, net, eqns, observed_CFs, idx_obj_map, α))
    (minf, minx, ret) = NLopt.optimize(opt, x0)
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
optimize_bls!(net::HybridNetwork, oCFs; kwargs...)::Float64 = optimize_bls!(net, find_quartet_equations(net)[1], oCFs; kwargs...)


"""
The objective function that is maximized during network optimization.
"""
function objective(X::Vector{T}, grad::Vector{T}, net::HybridNetwork, eqns::Array{QuartetData}, obsCFs::Matrix{T}, idx_obj_map::IdxObjMap, α::Float64)::T where T<:Float64
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
function setX!(net::HybridNetwork, X::Vector{Float64}, idx_obj_map::IdxObjMap)::Nothing
    for j in eachindex(X)
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


"""
Helper function to gather necessary information about `net` to
perform optimization.
"""
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
        if typeof(obj) <: Edge
            if getchild(obj).leaf continue end
            if getchild(obj).hybrid && getchild(getchild(obj)).leaf continue end
        end

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


"""
Helper function that takes a network `net` and its `param_map` (provided
by [`gather_optimization_info`](@ref)) and gathers each of the associated
parameters.
"""
function gather_params(net::HybridNetwork, param_map::Dict{Int, Int})::Array{Float64}
    params = zeros(length(param_map))
    for obj in vcat(net.hybrid, net.edge)
        if haskey(param_map, obj.number)
            params[param_map[obj.number]] = typeof(obj) <: Node ? getparentedgeminor(obj).gamma : obj.length
        end
    end
    return params
end


"""
Helper function that takes a network `net` and its `param_map` (provided
by [`gather_optimization_info`](@ref)) and gathers each of the associated
parameters.
"""
function gather_params(net::HybridNetwork)::Array{Float64}
    param_map = gather_optimization_info(net, true)[2]
    return gather_params(net, param_map)
end


"""
Computes the expected concordance factors of `net` with the inheritance
correlation parameter `α` (default=`Inf`).
"""
function compute_eCFs(net::HybridNetwork, α::Real=Inf)::Matrix{Float64}
    eqns, _, params, _ = find_quartet_equations(net)
    eCFs = zeros(length(eqns), 3)
    for j = 1:size(eCFs)[1]
        eCFs[j, 1], eCFs[j, 2] = compute_eCF(eqns[j], params, α)
        eCFs[j, 3] = 1 - eCFs[j, 1] - eCFs[j, 2]
    end
    return eCFs
end


"""
Creates a DataCF object containing the expected CFs for each quartet in `net`.
"""
function ExpectedDataCF(net::HybridNetwork, α::Real=Inf)::DataCF
    eqns, _, params, _ = find_quartet_equations(net);
    d = DataCF()
    for j in eachindex(eqns)
        eCF1, eCF2 = compute_eCF(eqns[j], params, α)
        push!(d.quartet, Quartet(j, eqns[j].q_taxa..., Vector{Float64}([eCF1, eCF2, 1.0 - eCF1 - eCF2])))
    end
    d.numQuartets = length(eqns)
    return d
end