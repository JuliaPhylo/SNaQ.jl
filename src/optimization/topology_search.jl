# Hill climbing over network topologies.

"""
    optimizetopology!(Nprime, old_eqns, move, params, qsub, q_idxs, opt_maxeval, force_resample_all, rng, ρ)

Optimizes the branch lengths and γ parameters of the network `Nprime`.

# Arguments:
- `Nprime::HybridNetwork`: network to be optimized, typically the proposed network in [`search`](@ref).
- `old_eqns::Vector{QuartetData}`: quartet equations of the network that came prior to `Nprime` in the
    hill climbing optimization of [`search`](@ref).
- `move::Symbol`: the topological move that turned the previous network into `Nprime`. Used to inform
    which quartets equations in `old_eqns` need to be re-calculated and which remain the same.
- `params::Tuple`: edges and/or hybrid nodes that define the move corresponding to `move` - used along with
    `move` for the aforementioned purpose.
- `qsub::AbstractMatrix{Float64}`: the observed quartet concordance factors of exactly the
    quartets in `q_idxs`, i.e. already sliced, with one row per entry of `old_eqns` and 3
    columns. Pre-sliced because `q_idxs` is fixed for a whole [`search`](@ref) run.
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
    qsub::AbstractMatrix{Float64},
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
        # `old_eqns[j]` and `Nprime_eqns[j]` are the same quartet, so its taxa are reused
        # rather than re-derived from `q_idxs[j]` (see `rebuildquartetequations!`).
        rebuildquartetequations!(Nprime, old_eqns, Nprime_eqns, Float64(ρ));
    end

    @debug "\tOptimizing branch lengths."
    Nprime_logPL = fitnumericalparameters!(Nprime, Nprime_eqns, qsub, ρ; maxeval=opt_maxeval, optargs...)

    return Nprime_logPL, Nprime_eqns
end


"""
    optimizetopology!(net, d)

This version is just a helper function for more clear tests. In the context of the
algorithm, this function recomputes values and wastes time.
"""
function optimizetopology!(net::HybridNetwork, d::DataCF)
    eqns = SNaQ.findquartetequations(net)[1];
    return fitnumericalparameters!(net, eqns, gatherCFmatrix(d); maxeval=500)
end


"""
    search(N, q, hmax; restrictions, ρ, propQuartets, propQuartetsFinal, preopt, probST, probQR, maxeval, maxequivPLs, opt_maxeval, seed, verbose, logfile)

Performs a single search for the optimal network topology with gradient-based optimization
of branch lengths and inheritance probabilities.

# Arguments
- `N::HybridNetwork`: The starting network topology.
- `q`: Observed quartet concordance factors.
- `hmax::Int`: Maximum number of hybridization events allowed.

# Optional Arguments
- `restrictions::Function=defaultrestrictions()`: Function to enforce restrictions on the proposed networks.
- `qinfTest::Bool`: whether to test for uninformative quartets (CFs near [1/3, 1/3, 1/3]) and
    exclude those quartets when performing the network search.
- `qtolAbs::Float64=1e-4`: the absolute tolerance used to detect uninformative quartets.
- `probQR::Float64=0.0`: probability at a given search iteration of utilizing weighted random
    sampling to (i) sample a poorly fitting quartet, (ii) sample an edge spanned by that
    quartet in the network, and finally (iii) only sample moves that include this edge.
- `ρ::Real=0.0`: inheritance correlation parameter in the range [0, 1]. `ρ = 0` corresponds to independent inheritance; `ρ = 1` corresponds to completely dependent inheritance.
- `propQuartets::Real=1.0`: Proportion of quartets to use during optimization.
- `preopt::Bool=false`: Whether to perform a pre-optimization step.
- `propQuartetsFinal::Real=1.0`: proportion of quartets in [0, 1] the final network is
  re-optimized on after the search. At 1.0 all quartets are used, but only when
  `propQuartets < 1.0` (otherwise the search already used them all). At 0.0 there is no
  re-optimization. The `propQuartetsFinal` sample used for re-optimization is drawn at
  random in each run, so [`multisearch`](@ref) re-scores the networks of all runs on one
  shared sample to fairly compare them.
- `probST::Real=0.3`: Probability of performing a subtree move before searching.
- `maxeval::Int=Int(1e8)`: Maximum number of evaluations.
- `maxequivPLs::Int=1500`: Maximum number of equivalent pseudo-likelihood scores to consider.
- `opt_maxeval::Int=max(30, N.numtaxa)`: Maximum evaluations when optimizing a proposed
  network's parameters. A proposal is fit from its parent's parameters, and the evaluations
  needed to converge that fit grow with the parameter count (`numtaxa - 3` for a tree), so a
  fixed budget under-converges large networks: at 200 taxa, 30 evaluations leave a residual
  3e4 times the acceptance gate and the proposal is rejected on its optimizer rather than on
  its topology.
- `seed::Int=abs(rand(Int) % 100000)`: Random seed for reproducibility.
- `verbose::Bool=false`: Whether to print verbose output.
- `logfile::String=""`: File to log detailed progress (used for debugging, but can also be used to examine convergence).

# Returns
- `best_network::HybridNetwork`: The network with the best (highest) composite log-likelihood.
- `best_score::Float64`: The composite log-likelihood of the best network.
"""
function search(
    N::HybridNetwork,
    q::Union{DataCF, AbstractMatrix{Float64}},
    hmax::Int;
    restrictions::Function=defaultrestrictions(),
    ρ::Real=0.0,
    propQuartets::Real=1.0,
    propQuartetsFinal::Real=1.0,
    preopt::Bool=true,
    probST::Real=0.3,
    probQR::Float64=0.0,
    maxeval::Int=Int(1e8),
    maxequivPLs::Int=1500,
    liktolAbs::Float64=1e-8,
    liktolRel::Float64=1e-4,
    opt_maxeval::Int=max(30, N.numtaxa),
    seed::Int=abs(rand(Int) % 100000),
    verbose::Bool=false,
    logfile::String="",
    filename::String="",
    outgroup::String="none",
    qinfTest::Bool=false,
    qtolAbs::Float64=1e-4,
    optargs...
)
    # Parameter enforcement
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")
    0 ≤ ρ ≤ 1 || error("ρ must be in range [0, 1] (ρ = $(ρ))")
    0 < propQuartets ≤ 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    0 ≤ propQuartetsFinal ≤ 1 || error("propQuartetsFinal must be in range [0, 1] (propQuartetsFinal = $(propQuartetsFinal))")
    0 ≤ probQR ≤ 1 || error("probQR must be in range [0, 1] (probQR = $(probQR))")
    0 ≤ probST ≤ 1 || error("probST must be in range [0, 1] (probST = $(probST))")
    outgroup == "none" || any(l -> l.name == outgroup, N.leaf) || error("No taxa in N have taxa name $(outgroup) (outgroup name)")
    qtolAbs ≥ 0.0 || error("qtolAbs must be ≥ 0.0 (qtolAbs = $qtolAbs)")
    0 ≤ probQR ≤ 1 || error("probQR must be in range [0, 1] (probQR = $probQR)")

    # Initial logging message
    starttime = time()
    filename != "" && open(string(filename, ".log"), "w+") do f end # clear any pre-existing text in the log file
    @logmessage filename """
    BEGIN: search with seed $(seed) at $(currenttime())
           starting topology: $(writenewick(N, round=true))"""

    # Convert q to a Matrix if it is a DataCF
    if typeof(q) <: DataCF
        q = gatherCFmatrix(q)
    end

    # Set the seed
    rng = Random.seed!(seed)

    # N = readnewick(writenewick(N));
    N = deepcopynetwork(N);
    for node in N.node
        if !node.leaf && !node.hybrid
            node.name = ""
        end
    end
    semidirectnetwork!(N)
    restrictions(N) || error("N does not meet restrictions IMMEDIATELY")

    if rand(rng) < probST
        found_different_net::Bool = false
        for j = 1:10_000
            try
                performrNNI1!(N, samplerNNIparameters(N, 1, rng)...);
                if restrictions(N)
                    found_different_net = true
                    break
                end
            catch
            finally
                if !found_different_net
                    # N = readnewick(writenewick(N));
                    N = deepcopynetwork(N);
                    semidirectnetwork!(N)
                end
            end
        end
        if !found_different_net
            @warn "Failed to adjust input network via NNI moves in a manner that met the provided restrictions (probST). Using the provided network instead."
        end
    end

    # Data used throughout the optimization process
    local q_idxs::Vector{Int64}
    if qinfTest
        informative = trues(size(q, 1))
        if qinfTest
            for (i, row) in enumerate(eachrow(q))
                informative[i] = isquartetinformative(row, qtolAbs)
            end
        end
        q_idxs = sampleqindices(N, propQuartets, informative, rng)
    elseif propQuartets == 1.0
        q_idxs = collect(1:nchoose4taxalength(N))
    else
        q_idxs = sampleqindices(N, propQuartets, rng)
    end
    current_logPL::Float64 = 0.0   # -logPL of `N`; the whole trace is never read, so it is not kept
    neq = findquartetequations(N, q_idxs);
    N_eqns::Vector{QuartetData} = neq[1];
    CFΔs::Vector{Float64} = []  # used when probQR != 0.0, computed WHEN NEEDED, so init'd to []
    unchanged_iters = 0

    # `q_idxs` is fixed for the whole run, so materialize its observed CFs once here
    # rather than re-slicing `q` every iteration -- costly for a lazy `q`.
    qsub::Matrix{Float64} = q[q_idxs, :]

    # Pre-optimizing the network's parameters
    if preopt
        @debug "Pre-optimizing"
        fitnumericalparameters!(N, N_eqns, qsub, ρ; optargs...)
        restrictions(N) || error("N does not meet restrictions after preopt")
        current_logPL = SNaQscore(N)
    else
        current_logPL = computeSNaQscore!(N_eqns, gatherparams(N), qsub, ρ)
    end

    moves_attempted = [];   # Vector of Tuples: (<move name>, <move parameters (i.e. nodes/edges)>)
    moves_proposed = Dict{Symbol,Int}()
    moves_accepted = Dict{Symbol,Int}()
    moves_logPL = Dict{Symbol,Vector{Float64}}()
    last_move = :nothing

    logtext(logfile, "Entering main loop with -logPL = $(current_logPL)")
    for j = 2:maxeval
        if j % 100 == 0
            logmoves(logfile, moves_proposed, moves_accepted, moves_logPL)
        end

        verbose && print("\rIteration $(j)/$(maxeval) - in a row=$(unchanged_iters)/$(maxequivPLs)              ")

        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        #Nprime = readnewick(writenewick(N));
        Nprime = deepcopynetwork(N);

        prop_move, prop_params = generatemoveproposal(Nprime, N_eqns, moves_attempted, hmax, probQR, qsub, CFΔs, rng, ρ)
        if prop_move === :none
            # Every move reachable from this topology has already been proposed and rejected
            logtext(logfile, "Iteration $(j): no untried move proposals remain, stopping.")
            break
        end
        last_move = prop_move
        applymove!(Nprime, prop_move, prop_params)
        @debug "Proposed move: $(prop_move), parameters: $(prop_params)"
        push!(moves_attempted, (prop_move, prop_params))
        @debug "Proposed: $(writenewick(Nprime, round=true))"

        if !haskey(moves_proposed, prop_move) moves_proposed[prop_move] = 0 end
        if !haskey(moves_accepted, prop_move) moves_accepted[prop_move] = 0 end
        if !haskey(moves_logPL, prop_move) moves_logPL[prop_move] = Vector{Float64}([]) end
        moves_proposed[prop_move] += 1

        # 2. Check for identifiability
        @debug "Proposed network level: $(getlevel(Nprime))"
        removedegree2nodes!(Nprime);
        while shrink3cycles!(Nprime) continue end
        while shrink2cycles!(Nprime) continue end   # keep shrinking until there is nothing to shrink

        # 2.2 Try re-rooting at the outgroup - if we can't, throw the network away
        if outgroup != "none"
            try
                rootatnode!(Nprime, outgroup)
            catch e
                if typeof(e) <: PN.RootMismatch
                    @debug "Nprime cannot be rooted at outgroup - skipping."
                    logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) (cannot reroot at outgroup)")
                    continue
                else
                    rethrow(e)
                end
            end
        end

        # 2.3 After removing some edges above, the root may have 2 edge now instead of 3 - we fix that here
        semidirectnetwork!(Nprime)

        # 3. Immediately throw away networks that don't meet restrictions 
        if !restrictions(Nprime)
            @debug "Nprime does not meet restrictions - skipping."
            logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) (restrictions not met)")
            continue
        end

        # Check whether we can do in-place updates here.
        # We CANNOT do inplace updates if:
        # 1. the number of hybrids changes, OR
        # 2. the number of optimization parameters in the network changed
        Nprime_np::Int = countoptimizationparams(Nprime)
        N_np::Int = countoptimizationparams(N)
        cannot_do_inplace::Bool = N.numhybrids != Nprime.numhybrids || N_np != Nprime_np

        # 4. Optimize branch lengths and compute logPL
        Nprime_logPL, Nprime_eqns = optimizetopology!(
            Nprime, N_eqns, prop_move, prop_params, qsub, q_idxs,
            opt_maxeval, cannot_do_inplace, rng, ρ; optargs...
        )
        Nprime_logPL == -Inf && error("Nprime_logPL is -Inf?? newick: $(writenewick(Nprime, round=true))\nold network: $(writenewick(N, round=true))\nprop move: $(prop_move)\nprop params: $(prop_params)")
        # computeSNaQscore!(Nprime, q) == Nprime_logPL || error("LOGPLS NOT EQUAL AFTER MOVE $(prop_move)")

        # 5. Accept / reject
        isnan(Nprime_logPL) && error("""
            Nprime_logPL = $(Nprime_logPL)
            $(writenewick(Nprime, round=true))
        """)
        if Nprime_logPL - current_logPL > liktolAbs && (current_logPL - Nprime_logPL) / current_logPL > liktolRel
            # Update current topology info
            N = Nprime
            N_eqns = Nprime_eqns
            CFΔs = []
            moves_accepted[prop_move] += 1
            push!(moves_logPL[prop_move], Nprime_logPL - current_logPL)
            current_logPL = Nprime_logPL

            # Log acceptance
            logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) ACCEPTED $(prop_move), new -logPL=$(round(current_logPL, digits=6))")

            # Update tracking vars
            unchanged_iters = 0
            moves_attempted = []
        else
            unchanged_iters += 1

            # Log rejection and reason
            logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) ($(round(Nprime_logPL, digits=3)) < $(round(current_logPL, digits=3)))")
        end

        # Early stopping checks
        if unchanged_iters > maxequivPLs
            @debug "stopping early after $(j) iterations"
            break
        end
    end
    SNaQscore!(N, current_logPL)

    if propQuartetsFinal == 1.0 && propQuartets != 1.0
        logmessage(filename, "Re-optimizing branch lengths with ALL quartets.")
        SNaQscore!(N, fitnumericalparameters!(N, q))
        logmessage(filename, "END propQuartets<1.0 post-search parameter optimization: found minimizer topology with SNaQ score=$(round(SNaQscore(N), digits=5))")
    elseif 0 < propQuartetsFinal < 1
        final_idxs = sampleqindices(N, propQuartetsFinal, rng)
        logmessage(filename, "Re-optimizing branch lengths with $(length(final_idxs)) of $(nchoose4taxalength(N)) quartets (propQuartetsFinal = $(propQuartetsFinal)).")
        final_eqns = findquartetequations(N, final_idxs)[1]
        SNaQscore!(N, fitnumericalparameters!(N, final_eqns, q[final_idxs, :], ρ))
        logmessage(filename, "END propQuartetsFinal<1.0 post-search parameter optimization: found minimizer topology with SNaQ score=$(round(SNaQscore(N), digits=5))")
    end

    # Remove internal node names that are not hybrids
    for node in N.node
        if !node.leaf && !node.hybrid
            node.name = ""
        end
    end

    # Rename hybrids to be from 1-H
    for (iH, H) in enumerate(N.hybrid)
        H.name = "H$iH"
    end

    logtext(logfile, "Search complete at $(currenttime()).\n\n")
    logmoves(logfile, moves_proposed, moves_accepted, moves_logPL)

    @logmessage filename "END: search with seed $(seed) after $(timeelapsed(time() - starttime)). SNaQ score = $(SNaQscore(N))"
    @logmessage filename writenewick(N)
    return N
end
