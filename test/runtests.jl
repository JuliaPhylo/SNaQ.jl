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
import SNaQ: writeExpCF, descData,
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
    numBad!, hasVeryBadTriangle!, index!, loglik!, blacklist!, cleaned!

import SNaQ: Node, findquartetequations, computeexpectedCFs,
    isvalidaddhybrid, computeexpectedDataCF,
    addhybrid!, removehybrid!,
    getparentedge, getparent,
    samplemovereticulateoriginparameters,
    isvalidmovereticulateorigin, movereticulateorigin!,deepcopynetwork,
    generatemoveproposal, applymove!,
    movereticulatetarget!, movereticulateorigin!,
    performrNNI!, performrNNI1!, performrNNI2!, performrNNI3!, performrNNI4!,
    isvalidrNNI1, isvalidrNNI2, isvalidrNNI3, isvalidrNNI4,
    allvalidrNNInodes, performrandomrNNI!,
    allvalidrNNI1nodes, allvalidrNNI2nodes,
    applymove!, performrSPR!, isvalidrSPR,
    samplerSPRparameters, multisearch, search

import PhyloNetworks: Node, Edge,
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
    printstyled(testfile, color=:blue)
    starttime = time()
    originalstdout = stdout
    originalstderr = stderr
    redirect_stdout(devnull)  # when testing for errors, comment these lines
    redirect_stderr(devnull)  # when testing for errors, comment these lines
    include(testfile)
    redirect_stdout(originalstdout)
    redirect_stderr(originalstderr)
    printstyled(": $(round(time() - starttime, digits=2))s elapsed\n", color=:black)
end


@testset "SNaQ.jl" begin
    runtestfile("test_bootstrap.jl")
    runtestfile("test_multipleAlleles.jl")
    runtestfile("test_perfectData.jl")
    runtestfile("test_readInputData.jl")

    runtestfile("test_propQuartets.jl")
    runtestfile("test_gradient_opt/test_opt_API.jl")
    runtestfile("test_gradient_opt/test_search_API.jl")
    runtestfile("test_gradient_opt/test_CF_recursive_blocks.jl")
    runtestfile("test_network_moves/test_add_remove_retic.jl")
    runtestfile("test_network_moves/test_move_target_origin.jl")
    runtestfile("test_network_moves/test_rNNI_moves.jl")
    runtestfile("test_network_moves/test_rSPR_moves.jl")
    runtestfile("test_network_moves/test_all_moves.jl")
end

