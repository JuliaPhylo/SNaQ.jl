using PhyloNetworks
include("../network_moves/add_remove_retic.jl")
include("../network_moves/rNNI_moves.jl")
include("../network_moves/rSPR_moves.jl")
include("../network_moves/move_origin_target.jl")


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
function search(N::HybridNetwork, q, hmax::Int; maxeval::Int=500, maxequivPLs::Int=100)
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")

    global API_WARNED
    if !API_WARNED
        @warn "API NOT YET FINALIZED"
        API_WARNED = true
    end

    N = deepcopy(N);
    if N.isrooted
        semidirect_network!(N)
        N.isrooted = false
    end

    logPLs = Array{Float64}(undef, maxeval)
    logPLs[1] = compute_logPL(N, q);
    prop_Ns = [N];
    unchanged_iters = 0

    for j = 2:maxeval
        print("\rCurrent best -logPL: $(-round(logPLs[j-1], digits=2))    ")

        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = readnewick(writenewick(N));
        propose_topology!(Nprime, hmax)
        @debug "Proposed: $(writenewick(Nprime, round=true))"

        # 2. Check for identifiability
        @debug "Proposed network level: $(getlevel(Nprime))"

        # 3. 2-cycles are NOT ALLOWED
        shrink2cycles!(Nprime);

        # 3. Optimize branch lengths
        @debug "\tgather quartets"
        q_eqns, _ = find_quartet_equations(Nprime)
        @debug "\toptimizing BLs"
        optimize_bls!(Nprime, q_eqns, q)
        @debug "\tdone optimizing BLs"

        max_gamma = -1.0
        for H in Nprime.hybrid
            max_gamma = max(max_gamma, max(getparentedge(H).gamma, getparentedgeminor(H).gamma))
        end
        min_BL = minimum(E.length for E in N.edge if !(E.hybrid && !E.ismajor))
        @debug "maxγ: $(round(max_gamma, digits = 2)) - minBL: $(round(min_BL, digits=4))"

        # 4. Compute logPL
        Nprime_logPL = compute_logPL(q_eqns, Nprime.edge, q)
        push!(prop_Ns, Nprime)

        # 5. Accept / reject
        if Nprime_logPL == NaN || abs(Nprime_logPL) < 1e-100
            error("Nprime_logPL = $(Nprime_logPL)")
        end

        if Nprime_logPL > logPLs[j-1]
            N = Nprime
            logPLs[j] = Nprime_logPL
            unchanged_iters = 0
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
    println("\rBest -logPL discovered: $(-round(logPLs[length(logPLs)], digits=2))")

    return N, logPLs

end



"""
Takes network `N` and modifies it with topological moves to generate a new proposal network.
"""
function propose_topology!(N::HybridNetwork, hmax::Int)

    # If any gammas are 0.0001 or 0.9999, propose removing that gamma
    bad_H = nothing
    for H in N.hybrid
        if getparentedgeminor(H).gamma <= 0.0001
            bad_H = H
            break
        end
    end

    if bad_H !== nothing
        @debug "Found H with γ=$(getparentedge(bad_H).gamma), removing it."
        remove_hybrid!(bad_H, N)
        return
    end

    if N.numhybrids < hmax && rand() < 0.50
        @debug "MOVE: add_random_hybrid!"
        add_random_hybrid!(N)
        return
    end

    if N.numhybrids > 0 && rand() < 0.50
        try
            if rand() < 0.5
                @debug "move_random_reticulate_origin!"
                move_random_reticulate_origin!(N)
                return
            else
                @debug "move_random_reticulate_target!"
                @debug writenewick(N, round=true)
                move_random_reticulate_target!(N)
                @debug writenewick(N, round=true)
                return
            end
        catch
            # Only catches if there were 0 valid moves for whichever above was chosen
        end
    end

    @debug "MOVE: perform_random_rNNI!"
    perform_random_rNNI!(N)
    return

end


