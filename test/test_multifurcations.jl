@testset "snaq! with a star tree as the starting topology" begin
    truetree = readnewick("((a:1,b:1):1,(c:1,d:1):1,(e:1,(f:1,g:1):1):1);")
    dcf = computeexpectedDataCF(truetree)
    star = readnewick("(a,b,c,d,e,f,g);")

    snaqnet = snaq!(star, dcf; hmax=0, Nfail=10, runs=10, seed=42, filename="")
    @test SNaQscore(snaqnet) > -1e-4
end

@testset "fitnumericalparameters! on a network with polytomies" begin
    net = readnewick("((((a:1,b:1,c:1):1)#H1:1::0.6,d:1):1,((#H1:1::0.4,e:1,f:1):1)#H2:1::0.7,((g:1)#H3:1::0.8,(h:1,#H3:1::0.2):1,#H2:1::0.3):1);")

    dcf = computeexpectedDataCF(net)
    Random.seed!(42)
    for E in net.edge
        E.length = rand() + 0.1
    end
    for H in net.hybrid
        γ = 0.5 + 0.4 * rand()
        getparentedge(H).gamma = γ
        getparentedgeminor(H).gamma = 1 - γ
    end

    polyfiterror = ErrorException("Found internal node in `N` with degrees 4. Resolve any polytomies manually (or with `resolvetreepolytomy!` from the `PhyloNetworks` package) before calling this function.")
    @test_throws polyfiterror fitnumericalparameters!(net, dcf; maxeval=10000)
end

@testset "snaq! with a star tree as the starting topology" begin
	truetree = readnewick("((a,b),(c,d),(e,(f,g)));");
	for E in truetree.edge E.length = 1.0 end
	dcf = computeexpectedDataCF(truetree);
	star = readnewick("(a,b,c,d,e,f,g);");
	snaqnet = snaq!(star, dcf; hmax=0, Nfail=20, runs=2, seed=42);
	@test snaqnet.numhybrids == 0
	@test SNaQscore(snaqnet) > -1e-4
end
