# Running independent searches and picking the best.

"""
    multisearch(
        N::HybridNetwork,
        q::Union{DataCF, AbstractArray{Float64}},
        hmax::Int;
        runs::Int=10,
        seed::Int=42,
        kwargs...
    ) -> Tuple{HybridNetwork, Vector{HybridNetwork}, Vector{Float64}}

Performs multiple independent searches to find the best network topology and parameters
that fit the observed quartet concordance factors.

# Arguments
- `N::Union{HybridNetwork, Vector{HybridNetwork}}`: The starting network topology or vector
    of topologies. If a vector, must be of length 1 or of length `runs`.
- `q::Union{DataCF, AbstractArray{Float64}}`: Observed quartet concordance factors.
- `hmax::Int`: Maximum number of hybridization events allowed.

# Optional Arguments
- `runs::Int=10`: Number of independent search runs.
- `seed::Int=42`: Random seed for reproducibility.
- `restrictions::Function=defaultrestrictions()`: Function that takes a `HybridNetwork` as
    its only argument and returns a `Bool`. Only networks that return `true` to this function
    will be considered during search.
- `propQuartetsFinal::Real=1.0`: passed to [`search`](@ref). If in (0, 1), each run's
    network is re-optimized on its own sample of quartets, so the networks of all runs are
    then re-scored on one shared sample (see [`rescoreonsharedquartets!`](@ref)).
- `ρ::Real=0.0`: inheritance correlation parameter, passed to [`search`](@ref).
- `kwargs...`: Additional keyword arguments passed to the [`search`](@ref) function.

# Returns
- `best_network::HybridNetwork`: The network with the best (highest) composite log-likelihood.
- `all_networks::Vector{HybridNetwork}`: All networks from the runs, sorted by score.
- `all_scores::Vector{Float64}`: Composite log-likelihood scores for each network.

# Notes
- This function uses distributed computing to perform searches in parallel.
- It returns the best network found across all runs.
"""
function multisearch(
    N::Union{HybridNetwork, AbstractVector{HybridNetwork}},
    q::Union{DataCF, AbstractMatrix{Float64}},
    hmax::Int;
    # Basic arguments
    verbose::Bool=true,
    runs::Int=10,
    seed::Int=abs(rand(Int) % 100000),
    logprefix::String="",
    filename::String="snaq",
    outgroup::String="none",
    restrictions::Function=defaultrestrictions(),
    propQuartetsFinal::Real=1.0,
    ρ::Real=0.0,
    kwargs...
)
    # Verify input parameters
    runs > 0 || error("runs must be > 0 (runs = $(runs)).")
    typeof(N) <: Vector{HybridNetwork} && length(N) != 1 && length(N) != runs && error("If N is a vector, it must be length 1 or length equal to runs (length(N) = $(length(N)), runs = $(runs))")
    N = deepcopy(N)

    if typeof(q) <: DataCF
            # If input is DataCF, make sure there's not a name mismatch
            dcf_names::Vector{String} = []
            for quartet in q.quartet
                for tax in quartet.taxon
                    if !(tax ∈ dcf_names)
                        push!(dcf_names, tax)
                    end
                end
            end

            # If any taxon are in the input but not the inputs but not the quartets, remove them from the inputs
            # If any taxon are NOT in the input but are in the quartets, error
            for T in (typeof(N) <: HybridNetwork ? [N] : N)
                any(tax -> tax ∉ tiplabels(T), dcf_names) && throw(ErrorException("Taxa in DataCF does not match taxa in input."))
                removetaxa = []
                for taxon in tiplabels(T)
                    if taxon ∉ dcf_names
                        push!(removetaxa, taxon)
                    end
                end
                if length(removetaxa) > 0
                    @warn "The following are in the inputs but not the DataCF, these taxa will be deleted: $(removetaxa)"
                    for taxon in removetaxa
                        PhyloNetworks.deleteleaf!(T, taxon)
                    end
                end
            end

            if typeof(N) <: HybridNetwork
                length(symdiff(tiplabels(N), dcf_names)) == 0 || throw(ErrorException("Taxa in DataCF does not match taxa in input."))
            else
                all(t -> length(symdiff(tiplabels(t), dcf_names)) == 0, N) || throw(ErrorException("Taxa in DataCF does not match taxa in one of the inputs."))
            end
        else
            # If input is AbstractArray{Float64}, make sure its size is (ntaxa choose 4, 3)
            nrow::Int = binomial(typeof(N) <: HybridNetwork ? N.numtaxa : N[1].numtaxa, 4)
            size(q) == (nrow, 3) || error("Input CF matrix should have size ($(nrow), 3), has size $(size(q)) instead.")
        end

    # Verify the starting network inputs
    Ns::Vector{HybridNetwork} = verifystartingtopologies!(N, outgroup, restrictions)

    # Convert q to a Matrix if it is a DataCF
    if typeof(q) <: DataCF
        q = gatherCFmatrix(q)
    end
    # Generate per-run seeds
    Random.seed!(seed)
    run_seeds = abs.(rand(Int, runs) .% 100000)

    filenames = begin
        # Log run details
        restrictionmsg = restrictions == defaultrestrictions() ? "default restrictions" :
            restrictions == restrictgallednetwork() ? "galled networks" :
            restrictions == restrictgalledtree() ? "galled trees" :
            restrictions == restrictrootedtreechild() ? "rooted tree child" :
            restrictions == restrictweaklytreechild() ? "weakly tree child" :
            restrictions == restrictstronglytreechild() ? "strongly tree child" :
            "custom restrictions"
        msg = """
        Beginning network optimization using SNaQ.jl across $runs runs with the following parameters:
            hmax = $hmax,
            seed = $seed,
            outgroup = $outgroup,
            restrictions = $restrictionmsg
        Root name for log files: $filename (absolute path $(abspath(filename)))
        Currently utilizing $(nprocs()) processor$(nprocs() > 1 ? "s" : "") and $(Threads.nthreads()) thread$(Threads.nthreads() > 1 ? "s" : "").
        """
        verbose && println(msg)

        if filename != ""
            open("$(filename).log", "w+") do f end
            logmessage(filename, msg)
        end

        # Make sure the log files can be created
        runs_path::String = ""
        if filename != ""
            runs_path = string(filename, "_runs/")
            mkpath(runs_path)
        end
        [filename == "" ? "" : "$(runs_path)run$(j)" for j = 1:runs]
    end

    # Do the runs distributed
    starttime = time()
    all_nets = pmap(
        j -> search(
            length(Ns) == 1 ? Ns[1] : Ns[j],
            q, hmax; seed = run_seeds[j], restrictions=restrictions,
            logfile = logprefix == "" ? "" : "$(logprefix)$(j)",
            filename = filenames[j],
            outgroup=outgroup, propQuartetsFinal=propQuartetsFinal, ρ=ρ, kwargs...
        ),
        1:runs
    )
    elapsed = timeelapsed(time() - starttime)

    if 0 < propQuartetsFinal < 1
        nused = rescoreonsharedquartets!(all_nets, q, propQuartetsFinal, seed, ρ)
        @logmessage filename "Re-scored the networks of all $runs runs on the same $nused quartets (propQuartetsFinal = $propQuartetsFinal)."
    end

    # Consolidate return data
    sort_idx = sortperm(SNaQscore.(all_nets), rev=true)
    bestnet = all_nets[sort_idx[1]]

    # Log results
    @logmessage filename """
    Finished optimizing topology at $(currenttime()) after $(elapsed).
    Optimal network: $(writenewick(bestnet, round=true))
    Optimal SNaQscore: $(SNaQscore(bestnet))
    To view all $runs inferred networks and their associated SNaQscore scores, see $(filename).out ($(abspath("$(filename).out")))"""

    if filename != ""
        open("$(filename).out", "w+") do f
            print(f,
                """
                $(writenewick(bestnet)) SNaQ score = $(SNaQscore(bestnet))
                Elapsed time: $(elapsed), $(runs) attempted runs

                -----------------------------------
                List of estimated networks for all runs (sorted by log-pseudolik; the larger, the better):
                """
            )
            for j in sort_idx
                println(f, " $(writenewick(all_nets[j])), with SNaQ score $(SNaQscore(all_nets[j]))")
            end
            println(f, "-----------------------------------")
        end

        open("$(filename).networks", "w+") do f
            for (j, i) in enumerate(sort_idx)
                write(f, "$(writenewick(all_nets[i])), with SNaQ score $(SNaQscore(all_nets[i]))")
                if j == 1
                    write(f, " (best network found, remaining sorted by log-pseudolik; the larger, the better)")
                end
                write(f, "\n")
            end
        end
    end

    # Clean up: the edges above roots have leftover values in them right now -
    #           we can't actually infer the lengths of those edges, so we clean
    #           those up here.
    for n in all_nets
        for L in n.leaf
            getparentedge(L).length = -1
        end
    end

    # Return
    return bestnet, all_nets[sort_idx]
end


"""
    rescoreonsharedquartets!(nets, q, propQuartetsFinal, seed, ρ) -> Int

Re-scores every network in `nets` on one sample of `propQuartetsFinal` of the quartets in `q`,
drawn with `seed`, without re-optimizing their parameters. Networks from independent runs are
optimized against different samples of quartets, so this makes their SNaQ scores comparable.
Returns the number of quartets in the shared sample.
"""
function rescoreonsharedquartets!(nets::Vector{HybridNetwork}, q::AbstractMatrix{Float64},
                                  propQuartetsFinal::Real, seed::Int, ρ::Real)::Int
    idxs = sampleqindices(size(q, 1), propQuartetsFinal, Random.seed!(seed))
    # Materialized once, not per network: every network is scored on the same sample.
    qsub::Matrix{Float64} = q[idxs, :]
    for net in nets
        eqns, _, params, _, _ = findquartetequations(net, idxs)
        SNaQscore!(net, computeSNaQscore!(eqns, params, qsub, Float64(ρ)))
    end
    return length(idxs)
end


"""
Verifies that starting topology(ies) `N` are ready to be optimized.
Modifies the network(s) `N` in-place.
"""
function verifystartingtopologies!(N::Union{HybridNetwork, AbstractVector{HybridNetwork}}, outgroup::String, restrictions::Function)::Vector{HybridNetwork}
    # Copy the input networks
    Ns::Vector{HybridNetwork} = typeof(N) <: HybridNetwork ? [deepcopynetwork(N)] : [deepcopynetwork(n) for n in N]
    for (j, n) in enumerate(Ns)
        # Split multifurcations
        if any(n -> length(n.edge) > 3, Ns[j].node)
            @warn "Input network #$(j) has a polytomy. SNaQ only infers binary networks, so this will be automatically resolved before inference."
            
            iters = 0
            while true
                iters += 1
                if iters > N.numnodes
                    error("Got stuck in an infinite loop while resolving polytomies. Please report this bug with your input tree(s) on GitHub.")
                end
                multi = findfirst(n -> length(n.edge) > 3, Ns[j].node)
                if isnothing(multi) break end
                PhyloNetworks.resolvetreepolytomy!(Ns[j], Ns[j].node[multi])
            end
        end

        # Prep data
        semidirectnetwork!(Ns[j]);

        # Make sure all leaf edges have some length so that code later doesn't error
        for E in Ns[j].edge
            E.length = E.length == -1.0 ? 0.0 : E.length
        end
        for H in Ns[j].hybrid
            if 1 ≥ getparentedge(H).gamma ≥ 0 && 1 ≥ getparentedgeminor(H).gamma ≥ 0 && getparentedge(H).gamma + getparentedgeminor(H).gamma ≈ 1
                continue
            end
            getparentedge(H).gamma = 0.5
            getparentedgeminor(H).gamma = 0.5
        end

        # Make sure starting network meets restrictions if any are provided
        restrictions(Ns[j]) || throw(ArgumentError("Starting topology #$(j) does not meet provided restrictions."))

        # If no outgroup exists, go next
        if outgroup == "none" continue end
        
        # Make sure the outgroup exists in this network
        if !any(L -> L.name == outgroup, Ns[j].leaf)
            throw(ArgumentError("Starting topology #$(j) does not contain the supplied outgroup ($(outgroup))."))
        end

        # Try rooting at outgroup if there is one
        try
            PN.rootatnode!(Ns[j], outgroup)
        catch e
            if typeof(e) <: PN.RootMismatch
                throw(ArgumentError("Starting topology #$(j) contains the outgroup but cannot be rooted at the outgroup."))
            else
                rethrow(e)
            end
        end
    end
    return Ns
end
