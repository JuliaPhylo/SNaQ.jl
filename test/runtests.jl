#if Threads.nthreads() == 1 error("Tests must be run with >1 thread.") end


using SNaQ
using PhyloNetworks
using Random
using DataFrames
using Distributed
using Test
using CSV
using Aqua
using StatsBase
using PhyloCoalSimulations
include("test_inplace_updates/misc.jl")

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

import SNaQ: Node, semidirect_network!, find_quartet_equations, compute_eCFs,
    is_valid_add_hybrid,
    add_hybrid!, remove_hybrid!,
    getparentedge, getparent,
    sample_move_reticulate_origin_parameters,
    is_valid_move_reticulate_origin, move_reticulate_origin!,deepcopy_network, semidirect_network!,
    generate_move_proposal, apply_move!,
    move_reticulate_target!, move_reticulate_origin!,
    semidirect_network!,
    perform_rNNI!, perform_rNNI1!, perform_rNNI2!, perform_rNNI3!, perform_rNNI4!,
    is_valid_rNNI1, is_valid_rNNI2, is_valid_rNNI3, is_valid_rNNI4,
    all_valid_rNNI_nodes, perform_random_rNNI!, semidirect_network!,
    all_valid_rNNI1_nodes, all_valid_rNNI2_nodes,
    apply_move!, perform_rSPR!, is_valid_rSPR, semidirect_network!,
    sample_rSPR_parameters

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
    # runtestfile("test_5taxon_readTopology.jl")
    # runtestfile("test_add2hyb.jl")
    # runtestfile("test_badDiamII.jl")
    # runtestfile("test_bootstrap.jl") 
    # runtestfile("test_calculateExpCF.jl")
    # runtestfile("test_calculateExpCF2.jl")
    # runtestfile("test_correctLik.jl")
    # runtestfile("test_deleteHybridizationUpdate.jl")
    # runtestfile("test_multipleAlleles.jl")
    # runtestfile("test_optBLparts.jl")
    # runtestfile("test_parameters.jl")
    # runtestfile("test_partition.jl")
    # runtestfile("test_partition2.jl")
    # runtestfile("test_perfectData.jl")
    # runtestfile("test_readInputData.jl")
    # runtestfile("test_hasEdge.jl")
    # runtestfile("test_undirectedOtherNetworks.jl")
    runtestfile("test_gradient_opt/test_opt_API.jl")
    runtestfile("test_gradient_opt/test_search_API.jl")
    runtestfile("test_gradient_opt/test_CF_recursive_blocks.jl")
    runtestfile("test_network_moves/test_add_remove_retic.jl")
    runtestfile("test_network_moves/test_move_target_origin.jl")
    runtestfile("test_network_moves/test_rNNI_moves.jl")
    runtestfile("test_network_moves/test_rSPR_moves.jl")
    runtestfile("test_network_moves/test_all_moves.jl")
end

