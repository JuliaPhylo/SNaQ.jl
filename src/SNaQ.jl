module SNaQ

    using Dates
    using Distributed

    using Printf: @printf
    using Random
    using Statistics: mean

    using Base.Threads

    # other libraries, indicate compatible version in Project.toml
    using CSV
    using DataFrames # innerjoin new in v0.21
    using DataStructures # for updateInCycle with priority queue
    using Distributions #for RateVariationAcrossSites
    using NLopt # for branch lengths optimization
    using StaticArrays
    using StatsBase # sample, etc.
    using PhyloNetworks

    import Base: show

    const DEBUGC = false # even more debug messages
    const PN = PhyloNetworks;
    const Edge = PN.Edge;
    const Node = PN.Node;
    global CHECKNET = false # for debugging only


    import PhyloNetworks: HybridNetwork, Network, Partition,
        tiplabels,
        isEqual, approxEq,
        assignhybridnames!, 
        getOtherNode, getIndex, getIndexEdge, getIndexNode,
        pushEdge!, pushNode!,
        setNode!, setEdge!,
        removeEdge!, removeNode!,
        addBL, deleteEdge!, deleteNode!,
        getconnectingedge, deleteIntNode!, numTreeEdges, numIntTreeEdges,
        hybridEdges, whichpartition, removeLeaf!, ##Almost used only in SNaQ functions
        pushHybrid!, removeHybrid!, printedges, printpartitions,
        samplebootstrap_multiloci, samplebootstrap_multiloci!, tree2Matrix,
        AQuartet, sort_stringasinteger!, tablequartetCF,
        RootMismatch

    export
        ## types & network definition
        DataCF,
        Quartet,
        # quartet CF
        readtrees2CF,
        readtableCF,
        readtableCF!,
        mapallelesCFtable,
        summarizedataCF,
        fittedquartetCF,
        # fitting: SNaQ and network bootstrap
        snaq!,
        readsnaqnetwork,
        bootsnaq,
        # functions to access relevant object variables
        loglik,
        loglik!,
        ########## New optimization functions
        optimize!,  # fixit: TODO: fix this name - it opts more than just BLs!
        computeloss,   # fixit: TODO: make this name in line with Julia conventions
        computeexpectedCFmatrix,  # fixit: TODO: no underscores in fxn names!!
        computeexpectedDataCF,             # fixit: TODO: again, not Julian naming convensions
        ########## New identifiability/restriction functions
        defaultrestrictions,    # fixit: TODO: some of these functions are called like
                                # restrictions=fxn(), while others are called like
                                # restirctions=fxn - standardize this
        norestrictions,
        knownidentifiable,
        restrictionset


    include("types.jl")
    include("auxiliary.jl")
    include("bootstrap.jl")
    include("multipleAlleles.jl")
    include("readquartetdata.jl")
    include("readwrite.jl")
    include("descriptive.jl")
    ############ NEW STUFF
    include("gradient_optimization/CF_recursive_blocks.jl")
    include("gradient_optimization/misc.jl")
    include("gradient_optimization/CF_equations.jl")
    include("gradient_optimization/inplace_updates.jl")
    include("network_properties/network_properties.jl")
    include("network_properties/identifiability_properties.jl")
    include("gradient_optimization/opt_API.jl")
    include("gradient_optimization/search_API.jl")
    include("gradient_optimization/BFS_search_distributed.jl")
    include("gradient_optimization/wrappers.jl")
    
    include("network_moves/misc.jl")
    include("network_moves/add_remove_retic.jl")
    include("network_moves/rNNI_validity.jl")
    include("network_moves/rNNI_moves.jl")
    include("network_moves/rSPR_validity.jl")
    include("network_moves/rSPR_moves.jl")
    include("network_moves/move_origin_target.jl")
    include("network_moves/identifiability_moves.jl")
    include("network_moves/flip_hybrid.jl")

end
