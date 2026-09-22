# Quickly generate a tree and its gene trees
Random.seed!(42)
tre = readnewick("(t1, t2, t3, t4, t5, t6, t7, t8);");
for E in tre.edge E.length = 0.0 end
tre = simulatecoalescent(tre, 1, 1)[1];
gts = simulatecoalescent(tre, 1000, 1);


@testset "LazyDataCF read from file can compute SNaQ score, but not with new CFs" begin
	ldcf = SNaQ.LazyDataCF(gts)
	s11 = computeSNaQscore!(tre, ldcf; propQuartets=0.1, seed=2)
	s12 = computeSNaQscore!(tre, ldcf; propQuartets=0.3, seed=12)
	s13 = computeSNaQscore!(tre, ldcf; numQuartets=15, seed=22)
	write("ldcf.csv", ldcf)

	ldcf = SNaQ.LazyDataCF("ldcf.csv")
	s21 = computeSNaQscore!(tre, ldcf; propQuartets=0.1, seed=2)
	s22 = computeSNaQscore!(tre, ldcf; propQuartets=0.3, seed=12)
	s23 = computeSNaQscore!(tre, ldcf; numQuartets=15, seed=22)

	@test s11 == s21
	@test s12 == s22
	@test s13 == s23
	rm("ldcf.csv")

	@test_throws Exception computeSNaQscore!(tre, ldcf; propQuartets=0.8, seed=0)
end

@testset "LazyDataCF can compute all SNaQ scores when it has all CFs" begin
	ldcf = SNaQ.LazyDataCF(gts)
	computeSNaQscore!(tre, ldcf; propQuartets=1.0)
	write("ldcf.csv", ldcf)

	ldcf = SNaQ.LazyDataCF("ldcf.csv")
	@test (computeSNaQscore!(tre, ldcf; numQuartets=40, seed=0); true)
	@test (computeSNaQscore!(tre, ldcf; numQuartets=70, seed=1); true)
	@test (computeSNaQscore!(tre, ldcf; propQuartets=0.25, seed=2); true)
	@test (computeSNaQscore!(tre, ldcf; propQuartets=1.0, seed=3); true)
	rm("ldcf.csv")
end
