module SNaQ

    using Dates
    using Distributed

    using Printf: @printf
    using Random
    using Statistics: mean

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
    global CHECKNET = false # for debugging only


    import PhyloNetworks: HybridNetwork, Edge, Node, Network, Partition,
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
        readnewicklevel1,
        readmultinewicklevel1,
        # quartet CF
        readtrees2CF,
        readtableCF,
        readtableCF!,
        readPhylip2CF,
        mapallelesCFtable,
        summarizedataCF,
        fittedquartetCF,
        # fitting: SNaQ and network bootstrap
        snaq!,
        readsnaqnetwork,
        topologymaxQpseudolik!,
        topologyQpseudolik!,
        bootsnaq,
        # functions to access relevant object variables
        loglik,
        loglik!,
        # Topological restrictions
        restrict_maximum_level,
        restrict_galled_tree,
        restrict_galled_network,
        restrict_rooted_tree_child,
        restrict_weakly_tree_child,
        restrict_strongly_tree_child,
        restriction_set


    include("types.jl")
    include("addHybrid_snaq.jl")
    include("auxiliary.jl")
    include("bootstrap.jl")
    include("deleteHybrid.jl")
    include("descriptive.jl")
    include("moves_snaq.jl")
    include("manipulateNet.jl")
    include("multipleAlleles.jl")
    include("parsimony.jl")
    include("pseudolik.jl")
    include("readquartetdata.jl")
    include("readwrite.jl")
    include("snaq_optimization.jl")
    include("undo.jl")
    include("update.jl")
    include("network_properties.jl")

end
