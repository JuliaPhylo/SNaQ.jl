using PhyloNetworks, SNaQ, DataFrames, PhyloCoalSimulations
using Test, Random
import SNaQ: multisearch

include(joinpath(@__DIR__, "../test_inplace_updates/misc.jl"))

### Network example
net = readnewick(joinpath(@__DIR__, "n1.netfile"))
q = SNaQ.computeexpectedCFs(net);
tre0 = majortree(net);


opt_rt = @elapsed opt_net, logPLs = multisearch(tre0, q, net.numhybrids; seed=7, maxequivPLs = 100)
@test computeloss(opt_net, q) > computeloss(tre0, q)

t0s = simulatecoalescent(tre0, 3, 1);
opt_rt = @elapsed opt_net, _ = multisearch(t0s, q, net.numhybrids; runs=length(t0s), seed=5, maxequivPLs = 100)
@test computeloss(opt_net, q) > minimum(computeloss(t0, q) for t0 in t0s)


@testset "Outgroups" begin
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
end

@testset "snaq! with propQuartets != 1.0" begin
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
	@test tcgidentifiable(snaqnet)
end

@testset "snaq! with custom restrictions" begin
	# snaq! with restrictions works as expected
	t0 = simulatecoalescent(tre0, 1, 1)[1];
	dcf = computeexpectedDataCF(net);
	snaqnet = snaq!(t0, dcf; restrictions=SNaQ.tcgidentifiable, hmax=4, runs=20, Nfail=50)
	@test SNaQ.tcgidentifiable(snaqnet)
end

# snaq! with restrictions works as expected
@testset "snaq! with restrictionset(;max_level=0) infers a tree, regardless of hmax" begin
	t0 = simulatecoalescent(tre0, 1, 1)[1];
	dcf = computeexpectedDataCF(net);
	snaqnet = snaq!(t0, dcf; restrictions=restrictionset(;max_level=0), hmax=4, runs=20, Nfail=50)
	@test snaqnet.numhybrids == 0
	@test tcgidentifiable(snaqnet)
end

@testset "snaq! respects restrictions when truth lies outside restrictions" begin
	# level-2 network topology
	truenet = readnewick("(((t4:0.4,(t8:0.3,t2:0.3):0.1):0.2,(t1:0.4,((t5:0.2,#H2:0.0::0.25):0.2,((t6:0.0,(t3:0.0)#H1:0.0::0.75):0.2)#H2:0.2::0.75):0.0):0.2):0.6,(t7:0.6,#H1:0.0::0.25):0.6);");
	Q = computeexpectedDataCF(truenet);
	T0 = simulatecoalescent(truenet, 1, 1)[1];
	snaqnet = snaq!(T0, Q; restrictions=restrictionset(;max_level=1), hmax=2, runs=3, Nfail=25)
	@test computeloss(snaqnet, Q) > computeloss(T0, Q)
	@test getlevel(snaqnet) == 1
end

@testset "snaq! with many hybrids does not error and improves the NCLL" begin
	for n in [10, 20]
		for h in [3, 6, 9]
			for seed in 1:10
				truenet = generate_net(n, h, seed)
				dcf = computeexpectedDataCF(truenet)
				T0 = simulatecoalescent(truenet, 1, 1)[1];
				rt = @elapsed snaqnet = snaq!(T0, dcf; restrictions=tcgidentifiable, propQuartets=0.1, hmax=h, runs=3, Nfail=5)
				@test computeloss(snaqnet, dcf) > computeloss(T0, dcf)
				@test tcgidentifiable(snaqnet)
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
	snaqnet = snaq!(net, dcf; hmax=net.numhybrids, Nfail=20, runs=3, opt_maxeval=10000, probST=0.0);
	@info computeloss(snaqnet, dcf)
	@test hardwiredclusterdistance(net, snaqnet, false) == 0
	@test computeloss(snaqnet, dcf) > -1e-2
end

@testset "snaq! with polytomies in the input" begin
	# Tree example
	multifurcation = readnewick("((a,b),(c,d),(e,f),(g,(h,i,j)));");
	for E in multifurcation.edge E.length = 1.0 end
	dcf = computeexpectedDataCF(multifurcation);
	for E in multifurcation.edge E.length = -1 end
	snaqnet = snaq!(multifurcation, dcf; hmax=0, Nfail=20, runs=3, probST=0.75);
	for E in reverse(snaqnet.edge)
		if !getchild(E).leaf && 0.0 ≤ E.length ≤ 1e-5	# hard-coded minimum for edge lengths
			PhyloNetworks.shrinkedge!(snaqnet, E)
		end
	end
	@test hardwiredclusterdistance(snaqnet, multifurcation, false) == 0

	# Network example
	multifurcation = generate_net(10, 3, 41);
	shrinkedges = sample([E for E in multifurcation.edge if !getchild(E).leaf && !E.hybrid], 5, replace=false)
	for E in shrinkedges
		PhyloNetworks.shrinkedge!(multifurcation, E)
	end
	dcf = computeexpectedDataCF(multifurcation)
	binary = SNaQ.verifystartingtopologies!(multifurcation, "none", (net) -> true)[1];
	computeloss(binary, dcf)
	snaqnet = snaq!(multifurcation, dcf; restrictions=norestrictions(), hmax=0, Nfail=20, runs=3, opt_maxeval=1000)
	@test loglik(snaqnet) > -1e-6
	for E in reverse(snaqnet.edge)
		if !getchild(E).leaf && !E.hybrid && 0.0 ≤ E.length ≤ 1e-5	# hard-coded minimum for edge lengths
			PhyloNetworks.shrinkedge!(snaqnet, E)
		end
	end
	@test hardwiredclusterdistance(snaqnet, multifurcation, false) == 0
end

@testset "snaq! with probQR = 0.5" begin
	tnet = generate_net(10, 2, 8);	# seed 8 is identifiable.
	dcf = computeexpectedDataCF(tnet);
	trueloss = computeloss(tnet, dcf); # approx 0, up to rounding error
	startnet = generate_net(10, 2, 12);
	startLL = computeloss(startnet, dcf)
	starthwcd = hardwiredclusterdistance(startnet, tnet, false)
	snaqnet = snaq!(startnet, dcf; hmax=tnet.numhybrids, Nfail=20, runs=3, opt_maxeval=100, probQR=0.5, probST=0.0, restrictions=SNaQ.norestrictions());
	@test computeloss(snaqnet, dcf) >= startLL
end

@testset "snaq! with probQR = 1.0" begin
	tnet = generate_net(10, 2, 8);	# seed 8 is identifiable.
	dcf = computeexpectedDataCF(tnet);
	trueloss = computeloss(tnet, dcf); # approx 0, up to rounding error
	startnet = generate_net(10, 2, 12);
	startLL = computeloss(startnet, dcf)
	starthwcd = hardwiredclusterdistance(startnet, tnet, false)
	snaqnet = snaq!(startnet, dcf; hmax=tnet.numhybrids, Nfail=20, runs=3, opt_maxeval=100, probQR=1.0, probST=0.0, restrictions=SNaQ.norestrictions());
	@test computeloss(snaqnet, dcf) >= startLL
end

@testset "snaq! with qinfTest=true" begin
	tnet = generate_net(10, 2, 8);	# seed 8 is identifiable.
	dcf = computeexpectedDataCF(tnet);
	trueloss = computeloss(tnet, dcf); # approx 0, up to rounding error
	startnet = generate_net(10, 2, 12);
	startLL = computeloss(startnet, dcf)
	starthwcd = hardwiredclusterdistance(startnet, tnet, false)
	snaqnet = snaq!(startnet, dcf; hmax=tnet.numhybrids, Nfail=20, runs=3, opt_maxeval=100, qinfTest=true, probST=0.0, restrictions=SNaQ.norestrictions());
	@test computeloss(snaqnet, dcf) >= startLL
end

@testset "snaq! with qinfTest=true, probQR=0.9, probST=1.0" begin
	tnet = generate_net(10, 2, 8);	# seed 8 is identifiable.
	dcf = computeexpectedDataCF(tnet);
	trueloss = computeloss(tnet, dcf); # approx 0, up to rounding error
	startnet = generate_net(10, 2, 12);
	startLL = computeloss(startnet, dcf)
	starthwcd = hardwiredclusterdistance(startnet, tnet, false)
	snaqnet = snaq!(startnet, dcf; hmax=tnet.numhybrids, Nfail=20, runs=3, opt_maxeval=100, qinfTest=true, probQR=0.9, probST=1.0, restrictions=SNaQ.norestrictions());
	@test computeloss(snaqnet, dcf) >= startLL
end

@testset "snaq! with ρ = 1.0" begin
	tnet = generate_net(10, 2, 8);	# seed 8 is identifiable.
	dcf = computeexpectedDataCF(tnet, 0.0);
	trueloss = computeloss(tnet, dcf, 0.0); # approx 0, up to rounding error
	startnet = generate_net(10, 2, 12);
	startLL = computeloss(startnet, dcf, 0.0)
	starthwcd = hardwiredclusterdistance(startnet, tnet, false)
	snaqnet = snaq!(startnet, dcf; hmax=tnet.numhybrids, Nfail=20, runs=3, opt_maxeval=100, ρ = 1.0, probST=0.0, restrictions=SNaQ.norestrictions());
	@test computeloss(snaqnet, dcf, 0.0) >= startLL
end

@testset "snaq! with ρ = 0.5" begin
	tnet = generate_net(10, 2, 8);	# seed 8 is identifiable.
	dcf = computeexpectedDataCF(tnet, 1.0);
	trueloss = computeloss(tnet, dcf, 1.0);
	startnet = generate_net(10, 2, 12);
	startLL = computeloss(startnet, dcf, 1.0)
	starthwcd = hardwiredclusterdistance(startnet, tnet, false)
	snaqnet = snaq!(startnet, dcf; hmax=tnet.numhybrids, Nfail=20, runs=3, opt_maxeval=100, ρ = 1.0, probST=0.0, restrictions=SNaQ.norestrictions());
	@test computeloss(snaqnet, dcf, 1.0) >= startLL
end
