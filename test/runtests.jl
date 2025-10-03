if Threads.nthreads() == 1 error("Tests must be run with >1 thread.") end


using SNaQ
using PhyloNetworks
using Random
using DataFrames
using Distributed
using Test
using CSV
using Aqua

## Import internal functions that are directly used in tests. There has got to be a better way
import SNaQ: checkNet,
    extractQuartet!, identifyQuartet!,
    eliminateHybridization!,
    updateSplit!, updateFormula!,
    calculateExpCF!, calculateExpCFAll!, logPseudoLik,
    updateInCycle!, updateContainRoot!, updateGammaz!,
    parameters!, update!, optTopRun1!,
    writeExpCF, writenewick_level1, descData,
    addHybridizationUpdate!, deleteHybridizationUpdate!,
    cleanBL!, cleanAfterRead!, identifyInCycle,
    updatePartition!, optBL!,undirectedOtherNetworks,hybridatnode!,
    # field getters
    istIdentifiable, fromBadDiamondI, inCycle, hasHybEdge,
    isBadDiamondI, isBadDiamondII, isExtBadTriangle, isVeryBadTriangle,
    k, typeHyb, gammaz,
    visited, edges_changed, nodes_changed, ht, numht,
    numBad, hasVeryBadTriangle, index, loglik, blacklist, cleaned,
    # field setters
    istIdentifiable!, fromBadDiamondI!, inCycle!, hasHybEdge!,
    isBadDiamondI!, isBadDiamondII!, isExtBadTriangle!, isVeryBadTriangle!,
    k!, typeHyb!, gammaz!,
    visited!, edges_changed!, nodes_changed!, ht!, numht!,
    numBad!, hasVeryBadTriangle!, index!, loglik!, blacklist!, cleaned!,
    chooseEdgesGamma, deleteHybrid!


import PhyloNetworks:
    Node, Edge,
    setNode!,setEdge!,
    approxEq,
    searchHybridNode,
    pushHybrid!

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(SNaQ;
    ambiguities = (broken = false),
    persistent_tasks = false,
    deps_compat = false)
end


include("test_output_helper.jl")
SNaQ.setCHECKNET(true)
@testset "SNaQ.jl" begin
    include("test_5taxon_readTopology.jl")
    include("test_add2hyb.jl")
    include("test_badDiamII.jl")
    include("test_bootstrap.jl") 
    include("test_calculateExpCF.jl")
    include("test_calculateExpCF2.jl")
    include("test_correctLik.jl")
    include("test_deleteHybridizationUpdate.jl")
    include("test_multipleAlleles.jl")
    include("test_optBLparts.jl")
    include("test_parameters.jl")
    include("test_partition.jl")
    include("test_partition2.jl")
    include("test_perfectData.jl")
    include("test_readInputData.jl")
    include("test_hasEdge.jl")
    include("test_undirectedOtherNetworks.jl")

    # Clean up some extra files that are created at some point in testing
    for filename in [
        "rand12Quartets.txt", "rand1Quartets.txt", "rand2Quartets.txt",
        "rand3Quartets.txt", "rand7Quartets.txt", "snaq.err", "snaq.log",
        "snaq.networks", "snaq.out", "summaryTreesQuartets.txt", "tableCF.txt"]
        if isfile(filename)
            rm(filename)
        end
    end
end
