using PhyloNetworks, SNaQ, DataFrames, PhyloCoalSimulations
using Test, Random

include(joinpath(@__DIR__, "../test_inplace_updates/misc.jl"))

### Network example
net = readnewick(joinpath(@__DIR__, "n1.netfile"))
q = SNaQ.computeexpectedCFs(net);
tre0 = majortree(net);


opt_rt = @elapsed opt_net, logPLs = multisearch(tre0, q, net.numhybrids; seed=7, maxequivPLs = 100)
@test hardwiredclusterdistance(net, opt_net, false) == 0

t0s = simulatecoalescent(tre0, 3, 1);
opt_rt = @elapsed opt_net, _ = multisearch(t0s, q, net.numhybrids; runs=length(t0s), seed=5, maxequivPLs = 100)
@test hardwiredclusterdistance(net, opt_net, false) == 0


# Test with outgroups
t0s = simulatecoalescent(tre0, 3, 1);
opt_rt = @elapsed opt_net, _ = multisearch(t0s, q, net.numhybrids; outgroup="a", runs=length(t0s), seed=5, maxequivPLs = 100)
@test any(t -> t.name == "a", getchildren(getroot(opt_net))) # NOT looking for 0 HWCD b/c "a" is not a true outgroup

t0s = simulatecoalescent(tre0, 3, 1);
opt_rt = @elapsed opt_net, _ = multisearch(t0s, q, net.numhybrids; outgroup="b", runs=length(t0s), seed=5, maxequivPLs = 100)
@test any(t -> t.name == "b", getchildren(getroot(opt_net))) # NOT looking for 0 HWCD b/c "b" is not a true outgroup

t0s = simulatecoalescent(tre0, 3, 1);
opt_rt = @elapsed opt_net, _ = multisearch(t0s, q, net.numhybrids; outgroup="e", runs=length(t0s), seed=5, maxequivPLs = 100)
@test hardwiredclusterdistance(net, opt_net, false) != 0 # "e" is NOT an outgroup, so we should have HWCD > 0
@test any(t -> t.name == "e", getchildren(getroot(opt_net)))

t0s = simulatecoalescent(tre0, 3, 1);
opt_rt = @elapsed opt_net, _ = multisearch(t0s, q, net.numhybrids; outgroup="f", runs=length(t0s), seed=5, maxequivPLs = 100)
@test hardwiredclusterdistance(net, opt_net, false) != 0 # "f" is NOT an outgroup, so we should have HWCD > 0
@test any(t -> t.name == "f", getchildren(getroot(opt_net)))

# Now test that the `snaq!` wrapper works as expected
t0 = simulatecoalescent(tre0, 1, 1)[1];
dcf = computeexpectedDataCF(net);
opt_rt = @elapsed snaqnet = snaq!(t0, dcf; hmax=1, propQuartets=0.85, runs=20, Nfail=50)
@test hardwiredclusterdistance(snaqnet, net, false) == 0
