# Quickly generate a tree and its gene trees
Random.seed!(42)
tre = readnewick("(t1, t2, t3, t4, t5, t6, t7, t8);");
for E in tre.edge E.length = 0.0 end
tre = simulatecoalescent(tre, 1, 1)[1];
gts = simulatecoalescent(tre, 1000, 1);


@testset "lazy DataCF read from file can compute SNaQ score, but not with new CFs" begin
	ldcf = DataCF(gts; lazy=true)
	s11 = computeSNaQscore!(tre, ldcf; propQuartets=0.1, seed=2)
	s12 = computeSNaQscore!(tre, ldcf; propQuartets=0.3, seed=12)
	s13 = computeSNaQscore!(tre, ldcf; numQuartets=15, seed=22)
	csvfile = tempname() * ".csv"
	write(csvfile, ldcf)

	ldcf = SNaQ.LazyDataCF(csvfile)
	@test ldcf isa DataCF && ldcf.lazy
	s21 = computeSNaQscore!(tre, ldcf; propQuartets=0.1, seed=2)
	s22 = computeSNaQscore!(tre, ldcf; propQuartets=0.3, seed=12)
	s23 = computeSNaQscore!(tre, ldcf; numQuartets=15, seed=22)

	@test s11 ≈ s21 atol=1e-10
	@test s12 ≈ s22 atol=1e-10
	@test s13 ≈ s23 atol=1e-10

	@test_throws ErrorException computeSNaQscore!(tre, ldcf; propQuartets=0.8, seed=0)
	@test_throws ErrorException ldcf.quartet[findfirst(i -> !haskey(ldcf.quartet, i), 1:70)]
	@test_throws ErrorException snaq!(tre, ldcf; hmax=0, propQuartets=0.5, propQuartetsFinal=0.5, runs=1, filename="")

	# with its gene trees, CFs missing from the file can be computed
	ldcf = SNaQ.LazyDataCF(csvfile, gts)
	@test computeSNaQscore!(tre, ldcf; propQuartets=0.8, seed=0) ≈
		computeSNaQscore!(tre, DataCF(gts; lazy=true); propQuartets=0.8, seed=0) atol=1e-10
	rm(csvfile)
end

@testset "lazy DataCF can compute all SNaQ scores when it has all CFs" begin
	ldcf = DataCF(gts; lazy=true)
	computeSNaQscore!(tre, ldcf; propQuartets=1.0)
	@test length(ldcf.quartet) == 70
	csvfile = tempname() * ".csv"
	write(csvfile, ldcf)

	ldcf = SNaQ.LazyDataCF(csvfile)
	@test (computeSNaQscore!(tre, ldcf; numQuartets=40, seed=0); true)
	@test (computeSNaQscore!(tre, ldcf; numQuartets=70, seed=1); true)
	@test (computeSNaQscore!(tre, ldcf; propQuartets=0.25, seed=2); true)
	@test (computeSNaQscore!(tre, ldcf; propQuartets=1.0, seed=3); true)
	rm(csvfile)
end

@testset "lazy DataCF is a LazyDataCF, and matches a DataCF with lazy=false" begin
	ldcf = DataCF(gts; lazy=true)
	@test ldcf.lazy && SNaQ.LazyDataCF(gts).lazy
	@test ldcf.numQuartets == 70
	@test length(ldcf.quartet) == 0
	@test tiplabels(ldcf) == ["t$i" for i in 1:8]
	@test computeSNaQscore!(tre, ldcf; numQuartets=20, seed=5) ≈
		computeSNaQscore!(tre, SNaQ.LazyDataCF(gts); numQuartets=20, seed=5) atol=1e-10
	@test length(ldcf.quartet) == 20
	@test [q.number for q in ldcf.quartet] == sort(collect(eachindex(ldcf.quartet)))

	# no propQuartets/numQuartets: scored on the CFs computed so far
	@test computeSNaQscore!(tre, ldcf) == SNaQscore(tre)
	@test length(ldcf.quartet) == 20

	@test computeSNaQscore!(tre, ldcf; propQuartets=1.0) ≈ computeSNaQscore!(tre, DataCF(gts)) atol=1e-10
	q = ldcf.quartet[1]
	@test q.taxon == ["t1", "t2", "t3", "t4"]
	@test q === ldcf.quartet[1]
	@test_throws BoundsError ldcf.quartet[71]

	# displaying a lazy DataCF or its quartets does not compute any CFs
	ldcf = DataCF(gts; lazy=true)
	sprint(show, MIME"text/plain"(), ldcf.quartet); sprint(show, ldcf.quartet); sprint(show, ldcf)
	@test length(ldcf.quartet) == 0
end

@testset "lazy DataCF safeguards" begin
	ldcf = DataCF(gts; lazy=true)
	@test_throws ErrorException DataCF(gts; lazy=true, whichQ="rand", numQ=10)
	@test_throws ErrorException computeSNaQscore!(tre, ldcf)	# no CFs computed yet
	@test_throws ErrorException computeSNaQscore!(tre, DataCF(gts); propQuartets=0.5)
	@test_throws ErrorException computeSNaQscore!(readnewick("(t1,t2,(t3,(t4,t5)));"), ldcf; propQuartets=0.5)
	@test_throws ErrorException write(tempname(), DataCF(gts))
	@test_throws ErrorException SNaQ.gatherCFmatrix(ldcf)
	@test_throws ErrorException fitnumericalparameters!(deepcopy(tre), ldcf)
	@test_throws ErrorException compositeloglik(tre, ldcf)
	@test_throws ErrorException compositedeviance(tre, ldcf)
	@test_throws ErrorException fittedquartetCF(ldcf)
	@test_throws ErrorException SNaQ.tablequartetCF(ldcf)
	@test_throws ErrorException summarizedataCF(ldcf)
	@test_throws ErrorException SNaQ.sorttaxa!(ldcf)
	@test_throws ErrorException SNaQ.taxadiff(ldcf, tre)
	@test_throws ErrorException multisearch(tre, ldcf, 0; runs=1, filename="")
	@test_throws ErrorException search(tre, ldcf, 0)
	@test length(ldcf.quartet) == 0	# none of the above computed any CF
end
