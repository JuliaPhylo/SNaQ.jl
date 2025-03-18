using PhyloNetworks, Random
include("../network_moves/add_remove_retic.jl")
include("../network_moves/rNNI_moves.jl")
include("../network_moves/rSPR_moves.jl")
include("../network_moves/move_origin_target.jl")
include("../network_moves/identifiability_moves.jl")
include("../network_properties/network_properties.jl")


API_WARNED = false
"""
API NOT FINALIZED
"""
function multi_search(N::HybridNetwork, q, hmax::Int; runs::Int=10, kwargs...)
    runs > 0 || error("runs must be > 0 (runs = $(runs)).")

    global API_WARNED
    if !API_WARNED
        @warn "API NOT YET FINALIZED"
        API_WARNED = true
    end

    N = readnewick(writenewick(N)); N.isrooted = false;
    perform_random_rNNI!(N, [1.0, 0.0, 0.0, 0.0])

    # Results from runs
    nets = Array{HybridNetwork}(undef, runs)
    logPLs = Array{Float64}(undef, runs)

    # Do the runs
    for j = 1:runs
        nets[j], all_logPLs = search(N, q, hmax; kwargs...)
        logPLs[j] = all_logPLs[length(all_logPLs)]
    end

    # Return
    sort_idx = sortperm(logPLs, rev=true)
    return nets[sort_idx[1]], nets[sort_idx], logPLs[sort_idx]
end


"""
API NOT FINALIZED
"""
function search(N::HybridNetwork, q, hmax::Int;
    restrictions::Function=no_restrictions(),
    α::Real=Inf,
    maxeval::Int=1500, maxequivPLs::Int=200)
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")

    global API_WARNED
    if !API_WARNED
        @warn "API NOT YET FINALIZED"
        API_WARNED = true
    end

    N = readnewick(writenewick(N));
    semidirect_network!(N)

    logPLs = Array{Float64}(undef, maxeval)
    N_qdata, _, N_params, _ = find_quartet_equations(N);
    logPLs[1] = compute_loss(N_qdata, N_params, q, α);
    prop_Ns = [N];
    unchanged_iters = 0

    moves_proposed = zeros(9)
    moves_accepted = zeros(9)

    for j = 2:maxeval
        # print("\rCurrent best -logPL: $(-round(logPLs[j-1], digits=2))          ($(j)/$(maxeval))                 ")

        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = readnewick(writenewick(N));
        move_proposed = propose_topology!(Nprime, hmax)
        moves_proposed[move_proposed] += 1
        @debug "Proposed: $(writenewick(Nprime, round=true))"

        # 2. Check for identifiability
        @debug "Proposed network level: $(getlevel(Nprime))"
        removedegree2nodes!(Nprime);
        # shrink3cycles!(Nprime);
        shrink2cycles!(Nprime);
        # shrink_bad_diamonds!(Nprime);

        # 3. Immediately throw away networks that don't meet restrictions
        if !restrictions(Nprime)
            @debug "Nprime does not meet restrictions - skipping."
            logPLs[j] = logPLs[j-1]
            continue
        end

        # 3. Optimize branch lengths
        @debug "\tgathering quartets"
        Nprime_qdata, _ = find_quartet_equations(Nprime);
        @debug "\toptimizing BLs"
        Nprime_logPL = optimize_bls!(Nprime, Nprime_qdata, q)
        @debug "Optimized network: $(writenewick(Nprime, round=true))"

        # 4. Remove hybrids with γ ≈ 0
        bad_H = nothing
        for H in N.hybrid
            if getparentedgeminor(H).gamma <= 0.001
                bad_H = H
                break
            end
        end
        if bad_H !== nothing
            @debug "Found H with γ=$(getparentedge(bad_H).gamma), removing it."
            remove_hybrid!(bad_H, N)
        end

        # 4. Compute logPL
        push!(prop_Ns, Nprime)

        # 5. Accept / reject
        if Nprime_logPL === NaN || abs(Nprime_logPL) < 1e-100
            error("Nprime_logPL = $(Nprime_logPL)")
        end

        if Nprime_logPL - logPLs[j-1] > 1e-8
            N = Nprime
            logPLs[j] = Nprime_logPL
            unchanged_iters = 0
            moves_accepted[move_proposed] += 1
        else
            logPLs[j] = logPLs[j-1]
            unchanged_iters += 1
        end

        # Early stopping checks
        if unchanged_iters > maxequivPLs
            @debug "stopping early after $(j) iterations"
            logPLs = logPLs[1:j]
            break
        end
    end
    # println("\rBest -logPL discovered: $(-round(logPLs[length(logPLs)], digits=2))")

    return N, logPLs

end


"""
Takes network `N` and modifies it with topological moves to generate a new proposal network.
Also updates `Nprime_eqns` in-place to hold the equations of the newly proposed network.
"""
function propose_topology!(N::HybridNetwork, hmax::Int)

    if N.numhybrids < hmax && rand() < 0.50
        @debug "MOVE: add_random_hybrid!"
        add_random_hybrid!(N)
        return 2
    end

    if rand() < 0.5

        @debug "MOVE: perform_random_rNNI!"
        return 2 + perform_random_rNNI!(N)
    
    else

        try
            if rand() < 0.5
                @debug "move_random_reticulate_origin!"
                # move_random_reticulate_origin!(N)
                move_random_reticulate_origin_local!(sample(N.hybrid), N)
                return 7
            else
                @debug "move_random_reticulate_target!"
                # move_random_reticulate_target!(N)
                move_random_reticulate_target_local!(sample(N.hybrid), N)
                return 8
            end
        catch
            # Only catches if there were 0 valid moves for whichever above was chosen
        end
    end

    @debug "MOVE: perform_random_rNNI!"
    return 2 + perform_random_rNNI!(N)

end


