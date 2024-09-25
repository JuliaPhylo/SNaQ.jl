using SNaQ
using PhyloNetworks
using Random
using DataFrames
using Test
using Aqua

@testset "SNaQ.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SNaQ;
        ambiguities = (broken=false),
        persistent_tasks=false)
    end
    
    SNaQ.setCHECKNET(true)

    @testset "SNaQ.jl" begin
        include("test_5taxon_readTopology.jl")
        include("test_add2hyb.jl")
        include("test_badDiamII.jl")


    end


end
