# PLACEHOLDER for now until all INCLUDE files are put into SNaQ.jl
using PhyloNetworks
for file in readdir(joinpath(@__DIR__, "../../src/network_moves/"), join=true)
    include(file)
end
include(joinpath(@__DIR__, "../../src/gradient_optimization/search_API.jl"))
include(joinpath(@__DIR__, "../test_inplace_updates/misc.jl"))