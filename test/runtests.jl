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
    deleteHybrid!, chooseEdgesGamma


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


function runtestfile(testfile::String)
    originalstdout = stdout
    originalstderr = stderr
    include(testfile)
    redirect_stdout(originalstdout)
    redirect_stderr(originalstderr)
end


SNaQ.setCHECKNET(true)
@testset "SNaQ.jl" begin
    runtestfile("test_5taxon_readTopology.jl")
    runtestfile("test_add2hyb.jl")
    runtestfile("test_badDiamII.jl")
    runtestfile("test_bootstrap.jl") 
    runtestfile("test_calculateExpCF.jl")
    runtestfile("test_calculateExpCF2.jl")
    runtestfile("test_correctLik.jl")
    runtestfile("test_deleteHybridizationUpdate.jl")
    runtestfile("test_multipleAlleles.jl")
    runtestfile("test_optBLparts.jl")
    runtestfile("test_parameters.jl")
    runtestfile("test_partition.jl")
    runtestfile("test_partition2.jl")
    runtestfile("test_perfectData.jl")
    runtestfile("test_readInputData.jl")
    runtestfile("test_hasEdge.jl")
    runtestfile("test_undirectedOtherNetworks.jl")
end

