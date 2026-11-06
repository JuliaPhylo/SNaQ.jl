# Proposing the next topology move.

"""
Checks whether a given quartet is informative based off of observed CFs.
Informative here is defined as any two entries in the quartet's observed
CFs having absolute difference greater than `atol`.
"""
function isquartetinformative(ocfrow::AbstractVector{Float64}, atol::Float64)
    for i = 1:2
        for j = (i+1):3
            if abs(ocfrow[i] - ocfrow[j]) > atol
                return true
            end
        end
    end
    return false
end


"""
Applies the move `move` on parameters `params` to network `N`.
"""
function applymove!(N::HybridNetwork, move::Symbol, params::Tuple)
    if move == :addhybrid
        return addhybrid!(N, params...)
    elseif move == :rNNI1
        return performrNNI1!(N, params...)
    elseif move == :rNNI2
        return performrNNI2!(N, params...)
    elseif move == :rNNI3
        return performrNNI3!(N, params...)
    elseif move == :rNNI4
        return performrNNI4!(N, params...)
    elseif move == :retic_origin || move == :retic_origin_local
        return movereticulateorigin!(N, params...)
    elseif move == :retic_target || move == :retic_target_local
        return movereticulatetarget!(N, params...)
    elseif move == :rSPR
        return performrSPR!(N, params...)
    elseif move == :flip_hybrid
        return fliphybrid!(N, params[1])
    end

    error("Move \"$(move)\" not recognized.")
end


"""
Randomly generates a move proposal from the function `samplemoveproposal`.
Sometimes the move sampled in `samplemoveproposal` will not have any valid
parameters, so this function repeatedly samples until a move with valid
parameters is selected. Also, makes sure the proposed move is not present in
`moves_attempted`, and appends the returned move to this vector.
"""
function generatemoveproposal(Nprime::HybridNetwork, N_eqns::Vector{QuartetData}, moves_attempted::Vector, hmax::Int, probQR::Float64, Q::Matrix{Float64}, CFΔs::Vector{Float64}, rng::TaskLocalRNG, ρ::Float64=0.0)::Tuple{Symbol,Any}
    required_edge::Union{Edge, Nothing} = nothing
    if probQR > 0.0 && rand(rng) <= probQR
        required_edge = sampleprobQRedge(Nprime, N_eqns, Q, CFΔs, rng, ρ)
    end

    validmove(mv::Symbol, pars) = !isnothing(pars) &&
        !alreadyattempted(moves_attempted, mv, pars) &&
        (isnothing(required_edge) || any(p -> p == required_edge, pars))

    retries::Int = 0
    move, params = samplemoveproposal(Nprime, hmax, rng)
    attempted_reqedges = 1

    while !validmove(move, params)
        move, params = samplemoveproposal(Nprime, hmax, rng)
        retries += 1
        if retries >= 1e3
            if isnothing(required_edge) || attempted_reqedges >= 100
                return (:none, ())      # caller stops the search; see `search`
            end
            # Sometimes we sample an edge that is not actually possible to
            # find in a move, so we re-sample another edge
            required_edge = sampleprobQRedge(Nprime, N_eqns, Q, CFΔs, rng, ρ)
            attempted_reqedges += 1
            retries = 0
        end
    end

    return (move, params)
end


generatemoveproposal(Nprime::HybridNetwork, ma::Vector, hmax::Int, rng::TaskLocalRNG) =
    generatemoveproposal(Nprime, Vector{QuartetData}([]), ma, hmax, 0.0, zeros(0, 0), zeros(0), rng, 0.0)


"""
Samples an edge from the network `N` with weights stored in `CFΔs`. If `CFΔs` is empty
(it starts empty and is reset to [] whenever the search finds a better network), these
weights are computed.
"""
function sampleprobQRedge(N::HybridNetwork, eqns::Vector{QuartetData}, Q::Matrix{Float64}, CFΔs::Vector{Float64}, rng::TaskLocalRNG, ρ::Float64=0.0)::Edge
    idxobjmap = gatheroptimizationinfo(N, false)[3]
    if length(CFΔs) == 0
        params = gatherparams(N);
        for (eqn, (ocf1, ocf2, ocf3)) in zip(eqns, eachrow(Q))
            ecf1, ecf2 = computeexpectedCF(eqn, params, ρ)
            ecf3 = 1.0 - ecf1 - ecf2
            push!(CFΔs, abs(ecf1 - ocf1) + abs(ecf2 - ocf2) + abs(ecf3 - ocf3))
        end
    end

    # Queue approach to run through every contributing equation in the sampled
    # CF's equation to make sure we select ALL edges that relate to this quartet
    iquartet = sample(rng, 1:length(eqns), Weights(CFΔs))
    Q = [eqns[iquartet].eqn]
    edges = []
    while length(Q) > 0
        curr = Q[1]
        deleteat!(Q, 1)

        append!(edges, curr.coal_edges)
        append!(Q, curr.divisions)
    end
    if length(edges) > 0
        return idxobjmap[sample(rng, unique(edges))]
    else
        validkeys = [k for k in keys(idxobjmap) if typeof(idxobjmap[k]) <: Edge]
        return idxobjmap[sample(validkeys)]
    end
end


"""
Helper function that determines whether the move `move` with parameters `params` has
already been attempted (i.e. is stored in the vector `moves_attempted`).
"""
function alreadyattempted(moves_attempted::Vector, move::Symbol, params::Tuple)::Bool
    for (amove, aparams) in moves_attempted
        move == amove || continue
        all_params_match::Bool = true
        for (aparam, param) in zip(aparams, params)
            if !(typeof(aparam) <: typeof(param))
                all_params_match = false
                break
            end
            if aparam.number != param.number || (typeof(aparam) <: Node && aparam.name != param.name)
                all_params_match = false
                break
            end
        end
        all_params_match && return true
    end
    return false
end


"""
Randomly samples a move to generate a new topology from `N`.
"""
function samplemoveproposal(N::HybridNetwork, hmax::Int, rng::TaskLocalRNG)::Tuple{Symbol,Any}
    if N.numhybrids < hmax && rand(rng) < 0.05
        @debug "SELECTED: add_random_hybrid!"
        return (:addhybrid, sampleaddhybridparameters(N, rng))
    end

    # If net has 0 hybrids, we can only do rNNI(1) or rSPR moves
    if N.numhybrids == 0
        # PROBABILITY OF EACH MOVE:
        # rNNI(1):  70%
        # rSPR:     30%

        r = rand(rng)
        if r < 0.7
            return (:rNNI1, samplerNNIparameters(N, 1, rng))
        else
            return (:rSPR, samplerSPRparameters(N, rng))
        end
    end


    # PROBABILITY OF EACH MOVE:
    # rNNI(1):      15%
    # rNNI(2):      0%  NEVER accepted!
    # rNNI(3):      0%  NEVER accepted!
    # rNNI(4):      10%
    # rSPR:         5%
    # origin:       10%
    # target:       10%
    # local origin: 20%
    # local target: 15%
    # fliphybrid:   15%
    probs = [0.15, 0.0, 0.0, 0.1, 0.05, 0.1, 0.1, 0.2, 0.15, 0.15]
    cumprobs = cumsum(probs)

    r = rand(rng)
    if r <= sum(cumprobs[1])
        return (:rNNI1, samplerNNIparameters(N, 1, rng))
    elseif r <= sum(cumprobs[2])
        return (:rNNI2, samplerNNIparameters(N, 2, rng))
    elseif r <= sum(cumprobs[3])
        return (:rNNI3, samplerNNIparameters(N, 3, rng))
    elseif r <= sum(cumprobs[4])
        return (:rNNI4, samplerNNIparameters(N, 4, rng))
    elseif r <= sum(cumprobs[5])
        return (:rSPR, samplerSPRparameters(N, rng))
    elseif r <= sum(cumprobs[6])
        return (:retic_origin, samplemovereticulateoriginparameters(N, rng))
    elseif r <= sum(cumprobs[7])
        return (:retic_target, samplemovereticulatetargetparameters(N, rng))
    elseif r <= sum(cumprobs[8])
        return (:retic_origin_local, samplemovereticulateoriginlocalparameters(N, rng))
    elseif r <= sum(cumprobs[9])
        return (:retic_target_local, samplemovereticulatetargetlocalparameters(N, rng))
    else
        return (:flip_hybrid, samplefliphybridparameters(N, rng))
    end

end



