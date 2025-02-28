using PhyloNetworks, Random
include("../network_moves/add_remove_retic.jl")
include("../network_moves/rNNI_moves.jl")
include("../network_moves/rSPR_moves.jl")
include("../network_moves/move_origin_target.jl")


function search(N::HybridNetwork, q, hmax::Int; maxeval::Int=500, maxequivPLs::Int=100, seed::Int=42)

    Random.seed!(seed)
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
        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = propose_topology(N, hmax)
        @debug "Proposed: $(writenewick(Nprime, round=true))"
        
        # 2. Check for identifiability


        # 3. Optimize branch lengths
        @info "\tgather quartets"
        q_eqns, _ = find_quartet_equations(Nprime)
        @info "\toptimizing BLs"
        optimize_bls!(Nprime, q_eqns, q)
        @info "\tdone optimizing BLs"

        max_gamma = -1.0
        for H in Nprime.hybrid
            max_gamma = max(max_gamma, max(getparentedge(H).gamma, getparentedgeminor(H).gamma))
        end
        min_BL = minimum(E.length for E in N.edge if !(E.hybrid && !E.ismajor))
        println("maxγ: $(round(max_gamma, digits = 2)) - minBL: $(round(min_BL, digits=4))")

        # 4. Compute logPL
        Nprime_logPL = compute_logPL(q_eqns, Nprime.edge, q)
        push!(prop_Ns, Nprime)

        # 5. Accept / reject
        if Nprime_logPL == NaN || abs(Nprime_logPL) < 1e-100
            error("Nprime_logPL = $(Nprime_logPL)")
        end

        if Nprime_logPL - logPLs[j-1] > 1e-4
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

    return N, logPLs

end



"""
Takes network `N` and modifies it with topological moves to generate a new proposal network.
"""
function propose_topology(N::HybridNetwork, hmax::Int)::HybridNetwork

    # If any gammas are 0.0001 or 0.9999, propose removing that gamma
    bad_H = nothing
    for H in N.hybrid
        if getparentedgeminor(H).gamma <= 0.0001
            bad_H = H
            break
        end
    end

    if bad_H !== nothing
        @info "Found H with γ=$(getparentedge(bad_H).gamma), removing it."
        Nprime = readnewick(writenewick(N)); Nprime.isrooted = false
        remove_hybrid!(bad_H, N)
        return Nprime
    end

    if N.numhybrids < hmax && rand() < 0.50
        @debug "MOVE: add_random_hybrid!"
        Nprime = readnewick(writenewick(N)); Nprime.isrooted = false
        add_random_hybrid!(Nprime)
        return Nprime
    end

    if N.numhybrids == 0 && rand() < 0.33

        @debug "MOVE: perform_random_rNNI!"
        Nprime = readnewick(writenewick(N)); Nprime.isrooted = false
        perform_random_rNNI!(Nprime)
        return Nprime
    
    else

        try
            if rand() < 0.5
                Nprime = readnewick(writenewick(N)); Nprime.isrooted = false
                @info "move_random_reticulate_origin!"
                move_random_reticulate_origin!(Nprime)
                return Nprime
            else
                Nprime = readnewick(writenewick(N)); Nprime.isrooted = false
                @info "move_random_reticulate_target!"
                @info writenewick(Nprime, round=true)
                move_random_reticulate_target!(Nprime)
                @info writenewick(Nprime, round=true)
                return Nprime
            end
        catch
            # Only catches if there were 0 valid moves for whichever above was chosen
        end
    end

    @debug "MOVE: perform_random_rNNI!"
    Nprime = readnewick(writenewick(N)); Nprime.isrooted = false
    perform_random_rNNI!(Nprime)
    return Nprime

end