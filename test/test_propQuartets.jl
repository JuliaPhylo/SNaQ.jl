# Tests to ensure that:
# 1. Final logliks are computed with ALL quartets, even with propQuartets < 1.0
# 2. Search finds improvements from initial topology when propQuartets < 1.0
#    (this was an old bug where propQuartets much smaller than 1.0 made it
#    impossible to improve the proposed topology during search)
# 3. A different set of quartets is chosen during each run when propQuartets < 1.0
# 4. propQuartetsFinal = 0.0 skips the final optimization, 0.0 < propQuartetsFinal < 1.0
#    re-optimizes each run on its own sample of quartets and then re-scores every run on
#    the same sample, and 1.0 is the default behavior
# 5. A lazy DataCF requires propQuartets < 1.0, and otherwise
#    finds the same networks as a DataCF with all observed CFs

truetre = readnewick("(((a,b):0.5,(c,d):0.5):0.5,(e,f):0.5);");
qcf = computeexpectedDataCF(truetre);
falsetre = readnewick("(((a,c):5.0,(e,b):5.0):1.0,(d,f):1.5);");


@testset "propQuartets<1.0: search improves over starting topology and final -loglik is computed with all quartets" begin
	snaqtre = snaq!(falsetre, qcf; hmax=0, propQuartets=0.5, seed=42, Nfail=10, runs=100, maxeval=5000);
	@test hardwiredclusterdistance(snaqtre, truetre, false) == 0
	@test loglik(snaqtre) ≈ computeSNaQscore!(truetre, qcf) atol=1e-10
end

@testset "propQuartets<1.0: different quartets sampled each run" begin
    rng1 = Random.seed!(1)
    idxs1 = SNaQ.sampleqindices(truetre, 0.5, rng1)
    rng2 = Random.seed!(2)
    idxs2 = SNaQ.sampleqindices(truetre, 0.5, rng2)
    @test idxs1 != idxs2
end

@testset "propQuartetsFinal=0.0: final optimization is skipped" begin
	snaqtre = snaq!(falsetre, qcf; hmax=0, propQuartets=0.5, propQuartetsFinal=0.0, seed=42, Nfail=10, runs=10, filename="pqf");
	@test hardwiredclusterdistance(snaqtre, truetre, false) == 0
	@test !occursin("Re-optimizing", read("pqf_runs/run1.log", String))
	rm("pqf_runs"; recursive=true); rm("pqf.log"); rm("pqf.out"); rm("pqf.networks")
end

@testset "0.0<propQuartetsFinal<1.0: every run is re-scored on the same sample of quartets" begin
	_, nets = multisearch(falsetre, qcf, 0; propQuartets=0.5, propQuartetsFinal=0.25, runs=3, seed=42, maxequivPLs=10, filename="pqf");
	idxs = SNaQ.sampleqindices(falsetre, 0.25, Random.seed!(42))
	Q = SNaQ.gatherCFmatrix(qcf)[idxs, :]
	for net in nets
		eqns, _, params, _, _ = findquartetequations(net, idxs)
		@test loglik(net) ≈ SNaQ.computeSNaQscore!(eqns, params, Q) atol=1e-10
	end
	@test occursin("propQuartetsFinal = 0.25", read("pqf_runs/run1.log", String))
	@test occursin("Re-scored the networks of all 3 runs", read("pqf.log", String))
	@test !occursin("Re-scored", read("pqf_runs/run1.log", String))
	rm("pqf_runs"; recursive=true); rm("pqf.log"); rm("pqf.out"); rm("pqf.networks")
end

@testset "propQuartetsFinal=1.0: same as the default" begin
	snaqtre1 = snaq!(falsetre, qcf; hmax=0, propQuartets=0.5, seed=42, Nfail=10, runs=5);
	snaqtre2 = snaq!(falsetre, qcf; hmax=0, propQuartets=0.5, propQuartetsFinal=1.0, seed=42, Nfail=10, runs=5);
	@test writenewick(snaqtre1) == writenewick(snaqtre2)
	@test loglik(snaqtre1) == loglik(snaqtre2)
end

@testset "propQuartetsFinal must be in [0, 1]" begin
	@test_throws Exception snaq!(falsetre, qcf; hmax=0, propQuartetsFinal=-0.1, runs=1)
	@test_throws Exception snaq!(falsetre, qcf; hmax=0, propQuartetsFinal=1.1, runs=1)
end

@testset "lazy DataCF" begin
	simtre = deepcopy(truetre)
	for E in simtre.edge if getchild(E).leaf E.length = 1.0 end end	# simulatecoalescent needs leaf edge lengths
	Random.seed!(42)
	gts = simulatecoalescent(simtre, 1000, 1);
	ldcf = DataCF(gts; lazy=true)
	@test_throws ErrorException snaq!(falsetre, ldcf; hmax=0)
	@test_throws ErrorException snaq!(falsetre, ldcf; hmax=0, propQuartets=1.0, propQuartetsFinal=0.5)
	@test_throws ErrorException computeSNaQscore!(truetre, ldcf)

	snaqtre = snaq!(falsetre, ldcf; hmax=0, propQuartets=0.8, propQuartetsFinal=1.0, seed=42, Nfail=10, runs=10, filename="");
	@test length(ldcf.quartet) == ldcf.numQuartets	# CFs computed by snaq! are kept in ldcf
	@test computeSNaQscore!(snaqtre, ldcf) ≈ loglik(snaqtre) atol=1e-10
	@test hardwiredclusterdistance(snaqtre, truetre, false) == 0
	@test loglik(snaqtre) ≈ computeSNaQscore!(snaqtre, gts) atol=1e-10
	@test computeSNaQscore!(truetre, DataCF(gts)) ≈ computeSNaQscore!(truetre, gts) atol=1e-10
end
