
const IdxObjMap = Dict{Int, Union{Node, Edge}};   # for readability


"""
    optimizetopology!(Nprime, old_eqns, move, params, q, q_idxs, opt_maxeval, force_resample_all, rng, ρ)

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
    I.e., when `propQuartets` in [`search`](@ref) is 1.0, this contains all integers from 1 to [`nchoose4taxalength`](@ref).
    When `propQuartets` is 0.1, it contains 10% as many integers, randomly selected in this range.
- `opt_maxeval::Int`: maximum number of evaluations when optimizing parameters.
- `force_resample_all::Bool`: if `true`, ignore `move` and `params` and re-calculate *every* quartet CF
    equation. Typically set to `true` when, e.g., a new reticulation is added to the network.
- `rng::TaskLocalRNG`: `TaskLocalRNG` object from which random numbers are generated. Ensures reproducibility.
- `ρ::Float64`: inheritance correlation parameter in the range [0, 1] used in calculating pseudo-likelihoods.
"""
function optimizetopology!(
    Nprime::HybridNetwork,
    old_eqns::Vector{QuartetData},
    move::Symbol,
    params::Tuple,
    q::Matrix{Float64},
    q_idxs::Vector{Int},
    opt_maxeval::Int,
    force_resample_all::Bool,
    rng::TaskLocalRNG,
    ρ::Float64=0.0;
    optargs...
)::Tuple{Float64, Vector{QuartetData}}
    Nprime_eqns::Vector{QuartetData} = Array{QuartetData}(undef, length(old_eqns))
    if !force_resample_all && can_update_inplace(move)
        @debug "\tGathering updated quartet equations."
        _, param_map, idxobjmap, _ = gatheroptimizationinfo(Nprime, true)
        updatequartetequations!(old_eqns, Nprime_eqns, Nprime, param_map, move, params, ρ)
    else
        @debug "\tGathering quartet equations."
        findquartetequations!(Nprime, q_idxs, Nprime_eqns);
    end

    @debug "\tOptimizing branch lengths."
    Nprime_logPL = optimize!(Nprime, Nprime_eqns, q[q_idxs,:], ρ; maxeval=opt_maxeval, optargs...)

    return Nprime_logPL, Nprime_eqns
end



"""
    optimize!(net, eqns, observed_CFs, ρ)

Optimizes the branch lengths of network `net` which is defined by quartet concordance factor
equations `eqns` based on observed quartet CFs `observed_CFs` under inheritance correlation
parameter `ρ`. `maxeval` adjusts the maximum number of loss evaluations in the optimization
process. This overloaded function is used as a helper method that can directly take
Vectors of `PhyloNetworks.QuartetT` as input so that the user doesn't need to handle
converting the input data.
"""
function optimize!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::AbstractVector{<:PhyloNetworks.QuartetT},
    ρ::Real=0.0;
    maxeval::Int=25
)::Float64
    obsCF_static = Array{Float64}(undef, length(observed_CFs), 3)
    for j in eachindex(observed_CFs)
        for k = 1:3
            obsCF_static[j, k] = observed_CFs[j].data[k]
        end
    end
    return optimize!(net, eqns, obsCF_static, ρ; maxeval=maxeval)
end

"""
Deprecated internal function - used for backwards compatibility in niche cases.
"""
optimize_bls!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::AbstractVector{<:PhyloNetworks.QuartetT},
    ρ::Real=0.0; kwargs...) = optimize!(net, eqns, observed_CFs, ρ; kwargs...)


"""
    optimizetopology!(net, d)

This version is just a helper function for more clear tests. In the context of the
algorithm, this function recomputes values and wastes time.
"""
function optimizetopology!(net::HybridNetwork, d::DataCF)
    eqns = SNaQ.findquartetequations(net)[1];
    return optimize!(net, eqns, gatherCFmatrix(d); maxeval=500)
end


"""
    optimize!(net, eqns, observed_CFs, ρ; maxeval, ftolRel, ftolAbs, xtolRel, xtolAbs)

Optimizes the branch lengths (and γ parameters) of network `net`.

# Required Arguments:
- `net::HybridNetwork`: the network to be optimized
- `eqns::Array{QuartetData}`: the set of quartet CF equations the define `net`
- `observed_CFs::Matrix{Float64}`: the corresponding set of observed CFs. Each row
    in this matrix MUST line up with each entry in `eqns`. If this is called through
    the [`search`](@ref) method, then this is handled.

# Optional Arguments:
- `ρ::Float64=0.0`: the inheritance correlation parameter in [0, 1]. `ρ = 0` is independent inheritance; `ρ = 1` is completely dependent. (Default=`0.0`)
- `maxeval::Int=25`: the maximum number of loss evaluations during optimization. (Default=`25`)
- `ftolRel::Float64=1e-8`: optimization parameter passed to the NLOpt.jl optimizer.
- `ftolAbs::Float64=1e-8`: optimization parameter passed to the NLOpt.jl optimizer.
- `xtolRel::Float64=1e-8`: optimization parameter passed to the NLOpt.jl optimizer.
- `xtolAbs::Float64=1e-8`: optimization parameter passed to the NLOpt.jl optimizer.
"""
function optimize!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::Matrix{Float64},
    ρ::Real=0.0;
    maxeval::Int=25,
    ftolRel::Float64=1e-12,
    ftolAbs::Float64=1e-12,
    xtolRel::Float64=1e-12,
    xtolAbs::Float64=1e-12
)::Float64

    0 ≤ ρ ≤ 1 || error("ρ must be in range [0, 1] (ρ = $ρ)")
    # Make sure there are no NaNs in the network's edge lengths
    # This is a bug that only seems to happen on Linux for some reason,
    # so it is hard for me to track down the source of the error
    for edge in net.edge
        if isnan(edge.length) edge.length = 0.0 end
        getchild(edge).leaf && continue
        edge.length = max(edge.length, 1e-5)    # starting optimization on a boundary can lead to failure
    end

    # Make sure γ values don't start on boundaries (this should never be
    # the case, but doing this check is free and may avoid errors)
    for H in net.hybrid
        γ = getparentedge(H).gamma
        γ = min(γ, 1.0 - 1e-5)
    end

    α = rhotoalpha(ρ)
    narg, param_map, idx_obj_map, params, LB, UB, init_steps = gatheroptimizationinfo(net, false)
    #opt = Opt(NLopt.LD_TNEWTON_PRECOND, narg)  # more accurate, but takes longer
    opt = Opt(NLopt.LD_LBFGS, narg)     # faster, but less accurate
    # OGNET = SNaQ.deepcopynetwork(net); # used in debugging

    opt.maxeval = maxeval
    opt.ftol_rel = ftolRel
    opt.ftol_abs = ftolAbs
    opt.xtol_rel = xtolRel
    opt.xtol_abs = xtolAbs

    initial_step!(opt, init_steps)
    opt.lower_bounds = LB
    opt.upper_bounds = UB

    x0::Vector{Float64} = [min(ub, val) for (ub, val) in zip(UB, params)]
    x0 = min.(UB .- 1e-12, x0)
    x0 = max.(LB .+ 1e-12, x0)
    NLopt.max_objective!(opt, (x, grad) -> objective(x, grad, net, eqns, observed_CFs, idx_obj_map, α))
    (maxf, maxx, ret) = NLopt.optimize(opt, x0)
    
    setX!(net, maxx, idx_obj_map)
    if maxf == -Inf
        error("Optimization error: maxf == -Inf")
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

    loglik!(net, maxf)
    return maxf
end
optimize!(net::HybridNetwork, oCFs; kwargs...)::Float64 = optimize!(net, findquartetequations(net)[1], oCFs; kwargs...)


"""
    optimize!(net::HybridNetwork, dcf::DataCF)

Optimizes the parameters of `net` with the quartet concordance factor data
in `dcf`. Returns the estimated likelihood of the network, which can also
be accessed later with `loglik(net)`.

### Parameters
- `ρ` is the inheritance correlation parameter in the range [0, 1] (default 0).
  `ρ = 0` corresponds to independent inheritance; `ρ = 1` corresponds to completely
  dependent inheritance.
- `maxeval` specifies the maximum number of optimization evaluations that the `NLopt`
  optimizer will perform under the hood (default 100).
"""
function optimize!(net::HybridNetwork, dcf::DataCF, ρ::Float64=0.0; maxeval::Int=100, kwargs...)::Float64
    0 ≤ ρ ≤ 1 || error("ρ must be between 0 and 1.")
    α = rhotoalpha(ρ)
    semidirectnetwork!(net)
    for E in net.edge
        E.length = max(E.length, 0.0)
    end
    for H in net.hybrid
        if getparentedge(H).gamma == -1 || getparentedgeminor(H).gamma == -1
            getparentedge(H).gamma = 0.5
            getparentedgeminor(H).gamma = 0.5
        end
    end
    eqns, parammap, parameters, _ = findquartetequations(net);
    obsCFs = gatherCFmatrix(dcf)
    optimize!(net, eqns, obsCFs, ρ; maxeval=maxeval, kwargs...)

    for q in dcf.quartet
        eqn = findquartetequations4taxa(net, q.taxon, parammap, ρ)
        expCF1, expCF2 = computeexpectedCF(eqn, parameters, ρ)
        q.expCF = [expCF1, expCF2, 1.0 - expCF1 - expCF2]
    end

    return loglik(net)
end


"""
WARNING: COMPLETELY EXPERIMENTAL AND UNSUPPORTED.
"""
function optimize_bls_staticγ!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::Matrix{Float64},
    ρ::Real=0.0;
    maxeval::Int=100,
    ftolRel::Float64=1e-12,
    ftolAbs::Float64=1e-12,
    xtolRel::Float64=1e-8,
    xtolAbs::Float64=1e-8
)::Float64

    α = rhotoalpha(ρ)
    narg, param_map, idx_obj_map, params, LB, UB, init_steps = gatheroptimizationinfo(net, false)

    # Figure out static parameters
    static_map::Dict{Int, Float64} = Dict(findfirst(idx -> idx_obj_map[idx] == hyb, 1:length(keys(idx_obj_map))) => getparentedgeminor(hyb).gamma for hyb in net.hybrid)
    @inline function parameterswithstatics(x::Vector{Float64})
        xout::Vector{Float64} = zeros(length(x) + length(static_map)) .- 1.0
        for idx in keys(static_map)
            xout[idx] = static_map[idx]
        end
        xidx = 1
        for idx in eachindex(xout)
            if xout[idx] == -1.0
                xout[idx] = x[xidx]
                xidx += 1
            end
        end
        any(xout .< 0.0) && error("Error mapping static parameters.")
        return xout
    end

    incl_idxs = sort([j for j in keys(idx_obj_map) if typeof(idx_obj_map[j]) <: Edge])
    narg_NLopt = narg - length(static_map)
    opt = Opt(NLopt.LD_LBFGS, narg_NLopt)     # faster, but less accurate
    opt.maxeval = maxeval
    opt.ftol_rel = ftolRel
    opt.ftol_abs = ftolAbs
    opt.xtol_rel = xtolRel
    opt.xtol_abs = xtolAbs

    initial_step!(opt, init_steps[incl_idxs])
    opt.lower_bounds = LB[incl_idxs]
    opt.upper_bounds = UB[incl_idxs]

    x0::Vector{Float64} = [min(ub / 2.0, val) for (ub, val) in zip(UB, params[incl_idxs])]
    x0 = min.(UB[incl_idxs] .- 1e-12, x0)
    x0 = max.(LB[incl_idxs] .+ 1e-12, x0)

    # Establish a map for fixing some static parameters

    NLopt.max_objective!(opt, (x, grad) -> objective_staticγ(parameterswithstatics(x), grad, net, eqns, observed_CFs, idx_obj_map, incl_idxs, α))

    (minf, minx, ret) = NLopt.optimize(opt, x0)
    if ret == :FAILURE
        @warn "ERROR: optimization returned :FAILURE"
    end
    setX!(net, parameterswithstatics(minx), idx_obj_map)

    return minf
end


"""
The objective function that is maximized during network optimization.
"""
function objective(X::Vector{T}, grad::Vector{T}, net::HybridNetwork, eqns::Array{QuartetData}, obsCFs::Matrix{T}, idx_obj_map::IdxObjMap, α::Float64)::T where T<:Float64
    setX!(net, X, idx_obj_map)
    fill!(grad, 0.0)
    loss = computelossandgradient!(eqns, X, grad, obsCFs, α)
    return loss
end


"""
WARNING: COMPLETELY EXPERIMENTAL AND UNSUPPORTED.

The objective function that is maximized during network optimization.
Provided γ value indices are static and therefore their gradients are not returned.
"""
function objective_staticγ(X::Vector{T}, grad::Vector{T}, net::HybridNetwork, eqns::Array{QuartetData}, obsCFs::Matrix{T}, idx_obj_map::IdxObjMap, include_idxs::Vector{Int}, α::Float64)::T where T<:Float64
    setX!(net, X, idx_obj_map)
    temp_grad = zeros(length(X))
    loss = computelossandgradient!(eqns, X, temp_grad, obsCFs, α)
    grad .= temp_grad[include_idxs]
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
function gatheroptimizationinfo(net::HybridNetwork, change_numbers::Bool=true)
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

            childnode = getchild(obj)
            if childnode.hybrid
                childnodechildren = getchildren(childnode)
                if length(childnodechildren) == 1 && childnodechildren[1].leaf
                    continue
                end
            end
            # if getchild(obj).hybrid && getchild(getchild(obj)).leaf continue end
        end

        haskey(param_map, obj.number) && error("Duplicate object number #$(obj.number).")
        param_map[obj.number] = param_idx
        idx_obj_map[param_idx] = obj
        param_idx += 1
    end

    params = gatherparams(net, param_map)
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
by [`gatheroptimizationinfo`](@ref)) and gathers each of the associated
parameters.
"""
function gatherparams(net::HybridNetwork, param_map::Dict{Int, Int})::Array{Float64}
    params = zeros(length(param_map))
    for obj in vcat(net.hybrid, net.edge)
        if haskey(param_map, obj.number)
            params[param_map[obj.number]] = typeof(obj) <: Node ? getparentedgeminor(obj).gamma : obj.length
        end
    end
    return params
end


"""
    gatherparams(net)

Helper function that takes a network `net` and its `param_map` (provided
by [`gatheroptimizationinfo`](@ref)) and gathers each of the associated
parameters.
"""
function gatherparams(net::HybridNetwork)::Array{Float64}
    param_map = gatheroptimizationinfo(net, true)[2]
    return gatherparams(net, param_map)
end


"""
    computeexpectedCFmatrix(net, ρ=0.0)

Computes the expected concordance factors of `net` with the inheritance
correlation parameter `ρ` (default=`0.0`). The returned `Matrix{Float64}`
object is unlabelled. See also [`computeexpectedDataCF`](@ref) for a `DataCF` object
with corresponding taxa information.
"""
function computeexpectedCFmatrix(net::HybridNetwork, ρ::Real=0.0)::Matrix{Float64}
    eqns, _, params, _ = findquartetequations(net)
    eCFs = zeros(length(eqns), 3)
    for j = 1:size(eCFs)[1]
        eCFs[j, 1], eCFs[j, 2] = computeexpectedCF(eqns[j], params, ρ)
        eCFs[j, 3] = 1 - eCFs[j, 1] - eCFs[j, 2]
    end
    return eCFs
end

"""
Deprecated - included for backwards compatibility in niche cases.
"""
computeexpectedCFs(net::HybridNetwork, ρ::Real=0.0)::Matrix{Float64} =
    computeexpectedCFmatrix(net, ρ)



"""
    computeexpectedDataCF(net, ρ=0.0)

Creates a DataCF object containing the expected CFs for each quartet in `net`.
"""
function computeexpectedDataCF(net::HybridNetwork, ρ::Real=0.0)::DataCF
    eqns, _, params, _ = findquartetequations(net);
    d = DataCF()
    for j in eachindex(eqns)
        eCF1, eCF2 = computeexpectedCF(eqns[j], params, ρ)
        q = Quartet(j, eqns[j].q_taxa..., Vector{Float64}([eCF1, eCF2, 1.0 - eCF1 - eCF2]))
        q.expCF = q.obsCF
        push!(d.quartet, q)
    end
    d.numQuartets = length(eqns)
    return d
end