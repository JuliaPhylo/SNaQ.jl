using PhyloNetworks, SNaQ, DataFrames, PhyloCoalSimulations
using Test, Random
import SNaQ: multisearch

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

# snaq! with restrictions works as expected
t0 = simulatecoalescent(tre0, 1, 1)[1];
dcf = computeexpectedDataCF(net);
snaqnet = snaq!(t0, dcf; restrictions=restrictionset(;max_level=1), hmax=1, propQuartets=0.85, runs=20, Nfail=50)
@test hardwiredclusterdistance(snaqnet, net, false) == 0
@test restrictionset(;max_level=1)(snaqnet)
@test knownidentifiable(snaqnet)

# snaq! with restrictions works as expected
t0 = simulatecoalescent(tre0, 1, 1)[1];
dcf = computeexpectedDataCF(net);
snaqnet = snaq!(t0, dcf; restrictions=SNaQ.knownidentifiable, hmax=4, runs=20, Nfail=50)
@test SNaQ.knownidentifiable(snaqnet)

# snaq! with restrictions works as expected
@testset "snaq! with restrictionset(;max_level=0) infers a tree, regardless of hmax" begin
	t0 = simulatecoalescent(tre0, 1, 1)[1];
	dcf = computeexpectedDataCF(net);
	snaqnet = snaq!(t0, dcf; restrictions=restrictionset(;max_level=0), hmax=4, runs=20, Nfail=50)
	@test snaqnet.numhybrids == 0
	@test knownidentifiable(snaqnet)
end

@testset "snaq! respects restrictions when truth lies outside restrictions" begin
	# level-2 network topology
	truenet = readnewick("(((t4:0.4,(t8:0.3,t2:0.3):0.1):0.2,(t1:0.4,((t5:0.2,#H2:0.0::0.25):0.2,((t6:0.0,(t3:0.0)#H1:0.0::0.75):0.2)#H2:0.2::0.75):0.0):0.2):0.6,(t7:0.6,#H1:0.0::0.25):0.6);");
	Q = computeexpectedDataCF(truenet);
	T0 = simulatecoalescent(truenet, 1, 1)[1];
	snaqnet = snaq!(T0, Q; restrictions=restrictionset(;max_level=1), hmax=2, runs=10, Nfail=25)
	@test computeloss(snaqnet, Q) > computeloss(T0, Q)
	@test getlevel(snaqnet) == 1
end

@testset "snaq! with larger networks" begin
	for n in [10, 15, 20]
		for h in [0, 3, 6, 9]
			for seed in 1:10
				truenet = generate_net(n, h, seed)
				dcf = computeexpectedDataCF(truenet)
				T0 = simulatecoalescent(truenet, 1, 1)[1];
				snaqnet = snaq!(T0, dcf; restrictions=SNaQ.knownidentifiable, propQuartets=0.25, hmax=h, runs=10, Nfail=50)
				@test computeloss(snaqnet, dcf) > computeloss(T0, dcf)
			end
		end
	end
end

@testset "snaq! works with missing edge lengths and gammas in starting topology" begin
	net = generate_net(10, 2, 8);	# seed 8 is identifiable.
	dcf = computeexpectedDataCF(net);
	trueloss = computeloss(net, dcf); # approx 0, up to rounding error
	for E in net.edge E.length = -1 end
	for H in net.hybrid
		getparentedge(H).gamma = -1
		getparentedgeminor(H).gamma = -1
	end
	snaqnet = snaq!(net, dcf; hmax=net.numhybrids, Nfail=20, runs=1, maxeval=10000, probST=0.0);
	@test hardwiredclusterdistance(net, snaqnet, false) == 0
	@test computeloss(snaqnet, dcf) ≈ trueloss atol=1e-12
end

@testset "snaq! with malformed inputs" begin
	multifurcation = readnewick("((a,b),(c,d),(e,f),(g,(h,i,j)));");
end
