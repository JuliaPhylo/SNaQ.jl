

"""
API NOT FINALIZED
"""
function multi_search(N::HybridNetwork, q, hmax::Int; runs::Int=10, seed::Int=abs(rand(Int) % 100000), kwargs...)
    # TODO: convert `q` to `DataCF` type and ensure that the taxa in that DataCF appear in the expected order
    runs > 0 || error("runs must be > 0 (runs = $(runs)).")

    # Prep data
    N = readnewick(writenewick(N));
    semidirect_network!(N)

    # Make sure starting network meets restrictions if any are provided
    (!haskey(kwargs, :restrictions) || kwargs[:restrictions](N)) || error("Starting topology does not meet provided restrictions.")

    # Generate per-run seeds
    Random.seed!(seed)
    run_seeds = abs.(rand(Int, runs) .% 100000)

    # Do the runs distributed
    nets_and_PLs = Distributed.pmap(1:runs) do j
        return search(N, q, hmax; seed = run_seeds[1], kwargs...)
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
    preopt::Bool=true,
    prehybprob::Real=0.3,
    prehybattempts::Int=5,
    probST::Real=0.3,
    maxeval::Int=Int(1e8),
    maxequivPLs::Int=1500,
    opt_maxeval::Int=10,
    seed::Int=abs(rand(Int) % 100000),
    verbose::Bool=false
)
    # Parameter enforcement
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")
    1 ≤ α ≤ Inf || error("α must be in range [1, ∞] (α = $(α))")
    0 < propQuartets ≤ 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    0 ≤ probST ≤ 1 || error("probST must be in range [0, 1] (probST = $(probST))")
    0 ≤ prehybprob ≤ 1 || error("prehybprob must be in range [0, 1] (prehybprob = $(prehybprob))")
    0 ≤ prehybattempts ≤ Inf || error("prehybattempts must be ≥ 0 (prehybattempts = $(prehybattempts))")

    # Set the seed
    rng = Random.seed!(seed)

    N = readnewick(writenewick(N));
    semidirect_network!(N)
    restrictions(N) || error("N does not meet restrictions IMMEDIATELY")

    if rand(rng) < probST
        try
            perform_rNNI1!(N, sample_rNNI_parameters(N, 1, rng)...);
        catch
        end
    end
    restrictions(N) || error("N does not meet restrictions after probST")

    # Initial pre-opt search
    if preopt
        @debug "Pre-optimizing"
        N, _ = search(N, q, hmax;
            restrictions=restrictions,
            preopt=false,
            probST=0.0,
            maxeval=1000,
            opt_maxeval=opt_maxeval,
            prehybattempts=0,
            maxequivPLs=50
        )
        restrictions(N) || error("N does not meet restrictions after preopt")
    end

    # Try many different ways of adding hybrids and pick the best
    if prehybattempts > 0 && rand(rng) < prehybprob
        @debug "Attempting pre-opt hybrid attachments"
        attempt_prehybs(N, q, hmax, restrictions, prehybattempts, rng)
        restrictions(N) || error("N does not meet restrictions after prehyb")
    end

    q_idxs = sample_qindices(N, propQuartets, rng)
    logPLs = Array{Float64}(undef, maxeval)
    logPLs[1], N_eqns = compute_loss(N, q, q_idxs, propQuartets, rng, α)
    unchanged_iters = 0

    moves_attempted = [];   # Vector of Tuples: (<move name>, <move parameters (i.e. nodes/edges)>)
    moves_proposed = Dict{Symbol,Int}()
    moves_accepted = Dict{Symbol,Int}()

    for j = 2:maxeval
        verbose && print("\rIteration $(j)/$(maxeval) - in a row=$(unchanged_iters)/$(maxequivPLs)              ")

        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = readnewick(writenewick(N));
        prop_move, prop_params = generate_move_proposal(Nprime, moves_attempted, hmax, rng)
        apply_move!(Nprime, prop_move, prop_params)
        push!(moves_attempted, (prop_move, prop_params))
        @debug "Proposed: $(writenewick(Nprime, round=true))"

        if !haskey(moves_proposed, prop_move) moves_proposed[prop_move] = 0 end
        if !haskey(moves_accepted, prop_move) moves_accepted[prop_move] = 0 end
        moves_proposed[prop_move] += 1

        # 2. Check for identifiability
        # TODO: check for galled-TC identifiability
        @debug "Proposed network level: $(getlevel(Nprime))"
        removedegree2nodes!(Nprime);
        while shrink3cycles!(Nprime) continue end
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
            moves_accepted[prop_move] += 1

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


function attempt_prehybs(net::HybridNetwork, q, hmax::Int, restrictions::Function, attempts::Int, rng::TaskLocalRNG)
    best_net = deepcopy_network(net)
    last_net = best_net
    best_PL = optimize_bls!(net, q)

    for nhyb = net.numhybrids:hmax
        best_PL = -Inf

        for j = 1:attempts
            local N0::HybridNetwork
            found_valid_move = false
            
            # Try to add a hybrid that complies w/ the given restrictions 100 times
            # If none are found, we just skip
            for restr_attempt = 1:100
                N0 = deepcopy_network(last_net)
                apply_move!(N0, :add_hybrid, Tuple(sample_add_hybrid_parameters(N0, rng)))
                if restrictions(N0)
                    found_valid_move = true
                    break
                end
            end
            found_valid_move || continue

            new_PL = optimize_bls!(N0, q)
            if new_PL > best_PL
                best_PL = new_PL
                best_net = N0
            end
        end

        best_net == last_net && break   # Couldn't find a move that is valid under `restrictions`
        last_net = best_net
    end
    return best_net, best_PL
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
    elseif move == :flip_hybrid
        return fliphybrid!(N, params[1])
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

    if N.numhybrids < hmax && rand(rng) < 0.15
        @debug "SELECTED: add_random_hybrid!"
        return (:add_hybrid, sample_add_hybrid_parameters(N, rng))
    end

    # If net has 0 hybrids, we can only do an NNI (SPRs are only
    # set up for reticulate moves at the moment)
    if N.numhybrids == 0
        return (:rNNI1, sample_rNNI_parameters(N, 1, rng))
    end


    # PROBABILITY OF EACH MOVE:
    # rNNI(1-4):  60%
    # moveretic:  20%
    # rSPR:       10%
    # fliphybrid: 10%
    r = rand(rng)
    if r < 0.60
        r = sample(rng, [1, 2, 3, 4])
        move_symbol = [:rNNI1, :rNNI2, :rNNI3, :rNNI4][r]
        return (move_symbol, sample_rNNI_parameters(N, r, rng))
    elseif r < 0.7
        return (:rSPR, sample_rSPR_parameters(N, rng))
    elseif r < 0.9
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
    else
        return (:flip_hybrid, sample_flip_hybrid_parameters(N, rng))
    end

end


