module SNaQ

    using Dates
    using Distributed
    using LinearAlgebra: diag, I, logdet, norm, LowerTriangular, mul!, lmul!, rmul!,
            Diagonal, cholesky, qr, BLAS
    # alternative: drop support for julia v1.4, because LinearAlgebra.rotate! requires julia v1.5
    # using LinearAlgebra # bring all of LinearAlgebra into scope
    # import LinearAlgebra.rotate! # allow re-definition of rotate!
    using Printf: @printf, @sprintf
    using Random
    using Statistics: mean, quantile, median

    using Base.Threads

    # other libraries, indicate compatible version in Project.toml
    using BioSequences
    using BioSymbols
    using Combinatorics: combinations
    using CSV
    using DataFrames # innerjoin new in v0.21
    using DataStructures # for updateInCycle with priority queue
    using Distributions #for RateVariationAcrossSites
    using FASTX
    using Functors: fmap
    using GLM # for the lm function
    using NLopt # for branch lengths optimization
    using StaticArrays
    using StatsBase # sample, coef etc.
    using StatsFuns # logsumexp, logaddexp, log2π, various cdf
    using StatsModels # re-exported by GLM. for ModelFrame ModelMatrix Formula etc
    using PhyloNetworks

    import Base: show
    import GLM: ftest
    import StatsModels: coefnames

    const DEBUGC = false # even more debug messages
    global CHECKNET = false # for debugging only


    import PhyloNetworks: HybridNetwork, tipLabels,getTipSubmatrix,
        resetNodeNumbers!,resetEdgeNumbers!,
        inheritanceWeight, isEqual,approxEq,
        assignhybridnames!, 
        getOtherNode,getIndex,getIndexEdge,getIndexNode,
        pushEdge!,pushNode!,
        setNode!,setEdge!,
        removeEdge!,removeNode!,
        searchHybridNode,searchHybridEdge,
        sampleBootstrapTrees, sampleBootstrapTrees!, tree2Matrix,
        addBL, deleteEdge!, deleteNode!,
        getconnectingedge, deleteIntNode!, numTreeEdges, numIntTreeEdges, ladderpartition,##. Possible PN jetsam. only used in SNaQ functions
        hybridEdges, whichPartition,removeLeaf!, ##Almost used only in SNaQ functions
        EdgeT, Node, Network, Partition


    export
        ## Network Definition
        DataCF,
        Quartet,
        readTopologyLevel1,
        tipLabels,
        sorttaxa!,
        ## SNAQ
        readTrees2CF,
        countquartetsintrees,
        readTableCF,
        readTableCF!,
        writeTableCF,
        readPhylip2CF,
        mapAllelesCFtable,
        readInputTrees,
        readnexus_treeblock,
        summarizeDataCF,
        snaq!,
        readSnaqNetwork,
        topologyMaxQPseudolik!,
        topologyQPseudolik!,
        ## Quartets
        Quartet,
        ## Network Manipulation
        fittedQuartetCF,
        ## Network Bootstrap
        bootsnaq,
        readMultiTopologyLevel1,
        ## parsimony
        parsimonySoftwired,
        parsimonyGF,
        maxParsimonyNet,
        readfastatodna
    


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
    include("readData.jl")
    include("readwrite.jl")
    include("snaq_optimization.jl")
    include("undo.jl")
    include("update.jl")

end
