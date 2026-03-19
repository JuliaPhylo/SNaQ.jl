# Tests to ensure that:
# 1. Final logliks are computed with ALL quartets, even with propQuartets < 1.0
# 2. Search finds improvements from initial topology when propQuartets < 1.0
#    (this was an old bug where propQuartets much smaller than 1.0 made it
#    impossible to improve the proposed topology during search)
# 3. A different set of quartets is chosen during each run when propQuartets < 1.0

truetre = readnewick("(((a,b):0.5,(c,d):0.5):0.5,(e,f):0.5);");
qcf = ExpectedDataCF(truetre);
falsetre = readnewick("(((a,c):5.0,(e,b):5.0):1.0,(d,f):1.5);");


@testset "propQuartets<1.0: search improves over starting topology and final -loglik is computed with all quartets" begin
	snaqtre = snaq!(falsetre, qcf; hmax=0, propQuartets=0.5, seed=42, Nfail=500, maxeval=5000);
	@test hardwiredclusterdistance(snaqtre, truetre, false) == 0
	@test loglik(snaqtre) ≈ compute_loss(truetre, qcf) atol=1e-10
end

@testset "propQuartets<1.0: different quartets sampled each run" begin
    rng1 = Random.seed!(1)
    idxs1 = SNaQ.sample_qindices(truetre, 0.5, rng1)
    rng2 = Random.seed!(2)
    idxs2 = SNaQ.sample_qindices(truetre, 0.5, rng2)
    @test idxs1 != idxs2
end
