# Fitting branch lengths and inheritance probabilities.

"""
    fitnumericalparameters!(net, eqns, observed_CFs, ρ)

Optimizes the branch lengths of network `net` which is defined by quartet concordance factor
equations `eqns` based on observed quartet CFs `observed_CFs` under inheritance correlation
parameter `ρ`. `maxeval` adjusts the maximum number of loss evaluations in the optimization
process. This overloaded function is used as a helper method that can directly take
Vectors of `PhyloNetworks.QuartetT` as input so that the user doesn't need to handle
converting the input data.
"""
function fitnumericalparameters!(
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
    return fitnumericalparameters!(net, eqns, obsCF_static, ρ; maxeval=maxeval)
end


"""
Deprecated internal function - used for backwards compatibility in niche cases.
"""
optimize_bls!(
    net::HybridNetwork,
    eqns::Array{QuartetData},
    observed_CFs::AbstractVector{<:PhyloNetworks.QuartetT},
    ρ::Real=0.0; kwargs...) = fitnumericalparameters!(net, eqns, observed_CFs, ρ; kwargs...)


"""
    fitnumericalparameters!(net, eqns, observed_CFs, ρ; maxeval, ftolRel, ftolAbs, xtolRel, xtolAbs)

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
function fitnumericalparameters!(
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
    # Built once per fit, then reused by all `maxeval` objective/gradient evaluations.
    batch = treequartetbatch(convert(Vector{QuartetData}, eqns), observed_CFs)
    NLopt.max_objective!(opt, (x, grad) -> objective(x, grad, net, batch, eqns, observed_CFs, idx_obj_map, α))
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

    SNaQscore!(net, maxf)
    return maxf
end


fitnumericalparameters!(net::HybridNetwork, oCFs; kwargs...)::Float64 = fitnumericalparameters!(net, findquartetequations(net)[1], oCFs; kwargs...)


"""
    fitnumericalparameters!(net, trees, ρ=0.0; propQuartets=1.0, seed=rand(Int), maxeval=100)
    fitnumericalparameters!(net, lazyq, ρ=0.0; propQuartets=1.0, seed=rand(Int), maxeval=100)

Optimizes the parameters of `net` directly against gene trees, or against a
[`LazyQuartetCF`](@ref) already built from them, without materializing the observed CFs of
every quartet. Returns the estimated likelihood, also readable with [`SNaQscore`](@ref).

`propQuartets` is the proportion of quartets to fit against, drawn with `seed`; the default
of 1.0 uses them all. Pass a smaller value for taxon counts where all `binomial(ntaxa,4)`
quartets will not fit.
"""
function fitnumericalparameters!(net::HybridNetwork, trees::Vector{HybridNetwork}, ρ::Real=0.0;
                                 propQuartets::Real=1.0, seed::Int=rand(Int), kwargs...)::Float64
    return fitnumericalparameters!(net, LazyQuartetCF(trees, sort(tiplabels(net))), ρ;
                                   propQuartets=propQuartets, seed=seed, kwargs...)
end

function fitnumericalparameters!(net::HybridNetwork, lazyq::LazyQuartetCF, ρ::Real=0.0;
                                 propQuartets::Real=1.0, seed::Int=rand(Int),
                                 maxeval::Int=100, kwargs...)::Float64
    0 ≤ ρ ≤ 1 || error("ρ must be between 0 and 1.")
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
    q_idxs, qsub = lazyquartetdata(net, lazyq, propQuartets, seed)
    eqns, _, _, _ = findquartetequations(net, q_idxs)
    return fitnumericalparameters!(net, eqns, qsub, ρ; maxeval=maxeval, kwargs...)
end


"""
    fitnumericalparameters!(net::HybridNetwork, dcf::DataCF)

Optimizes the parameters of `net` with the quartet concordance factor data
in `dcf`. Returns the estimated likelihood of the network, which can also
be accessed later with `SNaQscore(net)`.

### Parameters
- `ρ` is the inheritance correlation parameter in the range [0, 1] (default 0).
  `ρ = 0` corresponds to independent inheritance; `ρ = 1` corresponds to completely
  dependent inheritance.
- `maxeval` specifies the maximum number of optimization evaluations that the `NLopt`
  optimizer will perform under the hood (default 100).
"""
function fitnumericalparameters!(net::HybridNetwork, dcf::DataCF, ρ::Float64=0.0; maxeval::Int=100, kwargs...)::Float64
    0 ≤ ρ ≤ 1 || error("ρ must be between 0 and 1.")
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
    eqns, _, parameters, _ = findquartetequations(net);
    obsCFs = gatherCFmatrix(dcf)
    fitnumericalparameters!(net, eqns, obsCFs, ρ; maxeval=maxeval, kwargs...)

    for (i, q) in enumerate(dcf.quartet)
        eqn = eqns[i];
        expCF1, expCF2 = computeexpectedCF(eqn, parameters, ρ)
        q.expCF = [expCF1, expCF2, 1.0 - expCF1 - expCF2]
    end

    return SNaQscore(net)
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
[`objective`](@ref) using a prebuilt [`TreeQuartetBatch`](@ref).
"""
function objective(X::Vector{T}, grad::Vector{T}, net::HybridNetwork, batch::TreeQuartetBatch, eqns::Array{QuartetData}, obsCFs::Matrix{T}, idx_obj_map::IdxObjMap, α::Float64)::T where T<:Float64
    setX!(net, X, idx_obj_map)
    return computelossandgradient!(batch, convert(Vector{QuartetData}, eqns), X, grad, obsCFs, α)
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
