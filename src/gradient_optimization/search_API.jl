using PhyloNetworks, Random, Distributed
include("CF_equations.jl")
include("opt_API.jl")
include("../network_moves/add_remove_retic.jl")
include("../network_moves/rNNI_moves.jl")
include("../network_moves/rSPR_moves.jl")
include("../network_moves/move_origin_target.jl")
include("../network_moves/identifiability_moves.jl")
include("../network_properties/network_properties.jl")


"""
API NOT FINALIZED
"""
function multi_search(N::HybridNetwork, q, hmax::Int; runs::Int=10, seed::Int=abs(rand(Int) % 100000), kwargs...)
    runs > 0 || error("runs must be > 0 (runs = $(runs)).")

    # Prep data
    N = readnewick(writenewick(N));

    # Generate per-run seeds
    Random.seed!(seed)
    run_seeds = abs.(rand(Int, runs) .% 100000)

    # Do the runs distributed
    nets_and_PLs = Distributed.pmap(1:runs) do j
        return search(N, q, hmax; seed = run_seeds[j], kwargs...)
    end

    # Consolidate return data
    all_nets = Array{HybridNetwork}(undef, runs)
    all_logPLs = zeros(Float64, runs)
    for j = 1:runs
        all_nets[j], all_logPLs[j] = nets_and_PLs[j]
    end

    # Return
    sort_idx = sortperm(all_logPLs, rev=true)
    return all_nets[sort_idx[1]], all_nets[sort_idx], all_logPLs[sort_idx]
end


"""
API NOT FINALIZED
"""
function search(
    N::HybridNetwork, q, hmax::Int;
    restrictions::Function=no_restrictions(),
    α::Real=Inf,
    propQuartets::Real=1.0,
    probST::Real=0.3,
    maxeval::Int=Int(1e8),
    maxequivPLs::Int=200,
    opt_maxeval::Int=10,
    seed::Int=abs(rand(Int) % 100000)
)
    # Parameter enforcement
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")
    1 ≤ α ≤ Inf || error("α must be in range [1, ∞] (α = $(α))")
    0 < propQuartets ≤ 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    0 ≤ probST ≤ 1 || error("probST must be in range [0, 1] (probST = $(probST))")

    # Set the seed
    rng = Random.seed!(seed)

    N = readnewick(writenewick(N));
    semidirect_network!(N)

    if rand(rng) < probST
        try
            perform_rNNI1!(N, sample_rNNI_parameters(N, 1, rng)...);
        catch
        end
    end

    q_idxs = sample_qindices(N, propQuartets, rng)
    logPLs = Array{Float64}(undef, maxeval)
    logPLs[1], N_eqns = compute_loss(N, q, q_idxs, propQuartets, rng, α)
    unchanged_iters = 0

    moves_attempted = [];   # Vector of Tuples: (<move name>, <move parameters (i.e. nodes/edges)>)
    moves_proposed = zeros(10)
    moves_accepted = zeros(10)

    for j = 2:maxeval
        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = readnewick(writenewick(N));
        prop_move, prop_params = generate_move_proposal(Nprime, moves_attempted, hmax, rng)
        apply_move!(Nprime, prop_move, prop_params)
        push!(moves_attempted, (prop_move, prop_params))
        @debug "Proposed: $(writenewick(Nprime, round=true))"

        # 2. Check for identifiability
        # TODO: check for galled-TC identifiability
        @debug "Proposed network level: $(getlevel(Nprime))"
        removedegree2nodes!(Nprime);
        # shrink3cycles!(Nprime);
        while shrink2cycles!(Nprime) continue end   # keep shrinking until there is nothing to shrink
        # shrink_bad_diamonds!(Nprime);

        # 2.5 After removing some edges above, the root may have 2 edge now instead of 3 - we fix that here
        semidirect_network!(Nprime)

        # 3. Immediately throw away networks that don't meet restrictions
        if !restrictions(Nprime)
            @debug "Nprime does not meet restrictions - skipping."
            logPLs[j] = logPLs[j-1]
            continue
        end

        # 4. Optimize branch lengths and compute logPL
        Nprime_logPL, Nprime_eqns = optimize_topology!(
            Nprime, N_eqns, prop_move, prop_params, q, q_idxs,
            opt_maxeval, N.numhybrids != Nprime.numhybrids, rng, α
        )
        # compute_loss(Nprime, q) == Nprime_logPL || error("LOGPLS NOT EQUAL AFTER MOVE $(prop_move)")

        # 5. Accept / reject
        (Nprime_logPL === NaN || abs(Nprime_logPL) < eps()) && error("""
            Nprime_logPL = $(Nprime_logPL)
            $(writenewick(Nprime, round=true))
        """)
        if Nprime_logPL - logPLs[j-1] > 1e-8
            # Update current topology info
            N = Nprime
            N_eqns = Nprime_eqns
            logPLs[j] = Nprime_logPL

            # Update tracking vars
            unchanged_iters = 0
            moves_attempted = []
        else
            logPLs[j] = logPLs[j-1]
            unchanged_iters += 1
        end

        # 6. IF we chose Nprime (which is now N), remove hybrids with γ ≈ 0 from N
        #    if we did not choose Nprime, no point in doing this work
        if Nprime_logPL - logPLs[j-1] > 1e-8
            bad_Hs = []
            for H in N.hybrid
                if getparentedgeminor(H).gamma <= 0.001
                    bad_H = H
                end
            end
            if length(bad_Hs) != 0
                @error "Removing bad H!"
                @debug "Found $(length(bad_Hs)) hybrids with γ=$(getparentedge(bad_H).gamma), removing them."
                for H in bad_Hs
                    remove_hybrid!(N, H)
                end
                N_eqns = find_quartet_equations(N, q_idxs)
            end
        end

        # Early stopping checks
        if unchanged_iters > maxequivPLs
            @debug "stopping early after $(j) iterations"
            logPLs = logPLs[1:j]
            break
        end
    end

    return N, logPLs[length(logPLs)]

end


"""
Applies the move `move` on parameters `params` to network `N`.
"""
function apply_move!(N::HybridNetwork, move::Symbol, params::Tuple)
    if move == :add_hybrid
        return add_hybrid!(N, params...)
    elseif move == :rNNI1
        return perform_rNNI1!(N, params...)
    elseif move == :rNNI2
        return perform_rNNI2!(N, params...)
    elseif move == :rNNI3
        return perform_rNNI3!(N, params...)
    elseif move == :rNNI4
        return perform_rNNI4!(N, params...)
    elseif move == :retic_origin || move == :retic_origin_local
        return move_reticulate_origin!(N, params...)
    elseif move == :retic_target || move == :retic_target_local
        return move_reticulate_target!(N, params...)
    elseif move == :rSPR
        return perform_rSPR!(N, params...)
    end

    error("Move \"$(move)\" not recognized.")
end


"""
Randomly generates a move proposal from the function `sample_move_proposal`.
Sometimes the move sampled in `sample_move_proposal` will not have any valid
parameters, so this function repeatedly samples until a move with valid
parameters is selected. Also, makes sure the proposed move is not present in
`moves_attempted`, and appends the returned move to this vector.
"""
function generate_move_proposal(N::HybridNetwork, moves_attempted::Vector, hmax::Int, rng::TaskLocalRNG)
    retries::Int = 0
    move, params = sample_move_proposal(N, hmax, rng)

    while params === nothing || already_attempted(moves_attempted, move, params)
        move, params = sample_move_proposal(N, hmax, rng)
        retries += 1
        if retries >= 1e6
            error("Could not find any valid move proposals after 1e6 attempts.")
        end
    end

    return (move, params)
end


"""
Helper function that determines whether the move `move` with parameters `params` has
already been attempted (i.e. is stored in the vector `moves_attempted`).
"""
function already_attempted(moves_attempted::Vector, move::Symbol, params::Tuple)::Bool
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
        break
    end
    return false
end


"""
Randomly samples a move to generate a new topology from `N`.
"""
function sample_move_proposal(N::HybridNetwork, hmax::Int, rng::TaskLocalRNG)

    if N.numhybrids < hmax && rand(rng) < 0.05
        @debug "SELECTED: add_random_hybrid!"
        return (:add_hybrid, sample_add_hybrid_parameters(N, rng))
    end

    # If net has 0 hybrids, we can only do an NNI (SPRs are only
    # set up for reticulate moves at the moment)
    if N.numhybrids == 0
        return (:rNNI1, sample_rNNI_parameters(N, 1, rng))
    end

    r = rand(rng)
    if r < 0.66
        r = sample(rng, [1, 2, 3, 4])
        move_symbol = [:rNNI1, :rNNI2, :rNNI3, :rNNI4][r]
        return (move_symbol, sample_rNNI_parameters(N, r, rng))
    elseif r < 0.66
        return (:rSPR, sample_rSPR_parameters(N, rng))
    else
        r = sample(rng, [1, 2, 3, 4])   # 1 = move retic origin
                                        # 2 = move retic target
                                        # 3 = move retic origin local
                                        # 4 = move retic target local
        if r == 1
            return (:retic_origin, sample_move_reticulate_origin_parameters(N, rng))
        elseif r == 2
            return (:retic_target, sample_move_reticulate_target_parameters(N, rng))
        elseif r == 3
            return (:retic_origin_local, sample_move_reticulate_origin_local_parameters(N, rng))
        else
            return (:retic_target_local, sample_move_reticulate_target_local_parameters(N, rng))
        end
    end

end


