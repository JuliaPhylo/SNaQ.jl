using PhyloNetworks


function search(N::HybridNetwork, q, hmax::Int; maxeval::Int=500, maxequivPLs::Int=100)

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
        @info "Current: $(writenewick(N, round=true))"
        Nprime = propose_topology(N, hmax)
        @info "Proposed: $(writenewick(Nprime, round=true))"
        q_eqns, _ = find_quartet_equations(Nprime)
        optimize_bls!(Nprime, q_eqns, q)
        Nprime_logPL = compute_logPL(q_eqns, Nprime.edge, q)

        println("\t\t$(Nprime_logPL)")
        push!(prop_Ns, Nprime)

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
            @info "stopping early after $(j) iterations"
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

    if N.numhybrids < hmax && rand() < 0.20
        @info "MOVE: add_random_hybrid!"
        Nprime = deepcopy(N)
        add_random_hybrid!(Nprime)
        return Nprime
    end

    @info "MOVE: perform_random_rNNI!"
    Nprime = deepcopy(N)
    perform_random_rNNI!(Nprime)
    return Nprime

end