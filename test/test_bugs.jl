# Tests for bugs reported during PR #38:
#  1. `loglik` must equal the SNaQ v1 value (-100 * SNaQscore)
#  2. fix "custom restrictons" in the log files, when the restrictions were not custom
#  3. `compositedeviance` / `compositeloglik` / `q.expCF` fixes
#  4. `hmax` must be enforced on starting topologies with `N.numhybrids > hmax`
#  5. input processing can read v1 files and properly converts their scores to the new scale

bucky = readtableCF(joinpath(dirname(pathof(SNaQ)), "..", "examples", "buckyCF.csv"))
const NGENES_BUCKY = 30
truenewick = "((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);"
altnewick = "(D:0.0,C:0.0,(((E:0.0,#H1:9.773::0.215):0.563,O:0.0):0.515,(B:0.0,(A:9.773)#H1:0.0::0.785):18.178):19.157);"

scorefromexpCF(d::DataCF) = sum(q.obsCF[i] == 0 ? 0.0 : q.obsCF[i] * log(max(q.expCF[i], 1e-9) / q.obsCF[i])
                                for q in d.quartet for i in 1:3)

function cfbysplit(taxa::AbstractVector{<:AbstractString}, cf::AbstractVector{<:Real})
    d = Dict{Set{String},Float64}()
    for (k, j) in enumerate(2:4)
        side = Set{String}([taxa[1], taxa[j]])
        d[minimum(taxa) in side ? side : Set{String}(setdiff(taxa, side))] = cf[k]
    end
    return d
end

function expCFmatchesnetwork(d::DataCF, net::HybridNetwork; atol=1e-8)
    truth = Dict(Set(q.taxon) => cfbysplit(q.taxon, q.obsCF) for q in computeexpectedDataCF(net).quartet)
    for q in d.quartet
        stored = cfbysplit(q.taxon, q.expCF)
        expected = truth[Set(q.taxon)]
        all(isapprox(stored[k], expected[k]; atol=atol) for k in keys(expected)) || return false
    end
    return true
end

@testset "loglik keeps its v1 meaning, distinct from SNaQscore" begin
    net = readnewick(truenewick)
    s = computeSNaQscore!(net, bucky)
    @test SNaQscore(net) == s
    @test s < 0
    @test SNaQ.loglik(net) ≈ -100 * s
    @test SNaQ.loglik(net) > 0
    @test :loglik in names(SNaQ)
    SNaQ.loglik!(net, 50.0)
    @test SNaQscore(net) ≈ -0.5
    @test SNaQ.loglik(net) ≈ 50.0
end

@testset "expCF stored in DataCF is aligned with the network" begin
    d = deepcopy(bucky)
    net = readnewick(truenewick)
    s = computeSNaQscore!(net, d)
    @test scorefromexpCF(d) ≈ s
    @test expCFmatchesnetwork(d, net)

    s2 = fitnumericalparameters!(net, d)
    @test s2 > s
    @test scorefromexpCF(d) ≈ s2
    @test expCFmatchesnetwork(d, net)
    df = fittedquartetCF(d, :long)
    @test nrow(df) == 3 * length(d.quartet)

    d3 = deepcopy(bucky)
    est = snaq!(readnewick("((((D,C),(A,B)),E),O);"), d3; hmax=1, runs=1, Nfail=5, seed=123, filename="")
    @test scorefromexpCF(d3) ≈ computeSNaQscore!(est, deepcopy(bucky)) atol=1e-8
    @test expCFmatchesnetwork(d3, est)
end

@testset "compositedeviance and compositeloglik" begin
    d = deepcopy(bucky)
    obssum = sum(q.obsCF[i] == 0 ? 0.0 : q.ngenes * q.obsCF[i] * log(q.obsCF[i]) for q in d.quartet for i in 1:3)
    checkrelations(net, d, s) = begin
        dev = compositedeviance(net, d)
        @test dev ≈ -2 * NGENES_BUCKY * s
        @test compositeloglik(net, d) ≈ obssum - dev / 2
    end

    truenet = readnewick(truenewick)
    @test compositedeviance(truenet, deepcopy(bucky)) ≈ -2 * NGENES_BUCKY * computeSNaQscore!(truenet, deepcopy(bucky)) # fresh DataCF
    s = computeSNaQscore!(truenet, d)
    checkrelations(truenet, d, s)
    dev_before = compositedeviance(truenet, d)

    s2 = fitnumericalparameters!(truenet, d)
    checkrelations(truenet, d, s2)
    @test compositedeviance(truenet, d) < dev_before

    altnet = readnewick(altnewick)
    @test compositedeviance(altnet, d) ≈ -2 * NGENES_BUCKY * computeSNaQscore!(deepcopy(altnet), deepcopy(bucky))
    s3 = computeSNaQscore!(altnet, d)
    checkrelations(altnet, d, s3)
    s4 = fitnumericalparameters!(altnet, d)
    checkrelations(altnet, d, s4)

    obs = SNaQ.gatherCFmatrix(d)
    expmat = computeexpectedCFmatrix(altnet)
    @test compositeloglik(obs, expmat, NGENES_BUCKY) ≈ compositeloglik(altnet, d)
    @test compositeloglik(obs, expmat, fill(Float64(NGENES_BUCKY), size(obs, 1))) ≈ compositeloglik(altnet, d)
end

@testset "log file describes the restrictions used" begin
    tre = readnewick("((((D,C),(A,B)),E),O);")
    function restrictionline(; kw...)
        dir = mktempdir()
        fn = joinpath(dir, "snaq")
        snaq!(tre, bucky; hmax=0, runs=1, Nfail=1, seed=1, filename=fn, kw...)
        logtxt = read(fn * ".log", String)
        m = match(r"restrictions = (.*)", logtxt)
        return isnothing(m) ? "" : m.captures[1]
    end
    mycheck(net) = true
    for line in (restrictionline(), restrictionline(restrictions=norestrictions),
                 restrictionline(restrictions=norestrictions()),
                 restrictionline(restrictions=defaultrestrictions()))
        @test !occursin("custom", line)
        @test occursin("none", line)
    end
    for line in (restrictionline(restrictions=restrictgalledtree), restrictionline(restrictions=restrictgalledtree()))
        @test occursin("galled tree", line)
    end
    @test occursin("tree-child", restrictionline(restrictions=restrictstronglytreechild))
    @test occursin("level ≤ 2", restrictionline(restrictions=restrictmaximumlevel(2)))
    line = restrictionline(restrictions=restrictionset(max_level=2, galled_network=true))
    @test occursin("level ≤ 2", line) && occursin("galled network", line)
    @test occursin("mycheck", restrictionline(restrictions=mycheck))
end

@testset "hmax enforced on starting topologies" begin
    starth2 = readnewick("(((((D,C))#H2,((A,B))#H1),(#H1,E)),(#H2,O));")
    @test starth2.numhybrids == 2
    for hmax in (0, 1)
        est = @test_logs (:warn, r"hmax") match_mode=:any snaq!(starth2, bucky;
            hmax=hmax, runs=2, Nfail=5, seed=42, filename="")
        @test est.numhybrids <= hmax
    end
    @test starth2.numhybrids == 2
    starth2gam = readnewick("(((((D:1,C:1):1)#H2:1::0.8,((A:1,B:1):1)#H1:1::0.7):1,(#H1:1::0.3,E:1):1):1,(#H2:1::0.2,O:1):1);")
    @test_throws r"hmax" SNaQ.search(starth2gam, bucky, 0; filename="")
end

@testset "reading SNaQ v1 output files" begin
    v1net = readsnaqnetwork(joinpath(dirname(pathof(SNaQ)), "..", "examples", "net1_snaq.out"))
    @test SNaQ.loglik(v1net) ≈ 126.74250981974959
    @test SNaQscore(v1net) ≈ -1.2674250981974959
    v2net = readsnaqnetwork(joinpath(dirname(pathof(SNaQ)), "..", "examples", "net1.out"))
    @test SNaQscore(v2net) ≈ -0.0142092522636232
    fn = joinpath(mktempdir(), "v1.networks")
    write(fn, """
        (A,(B,(C,D)));, with -loglik 12.5 (best network found, remaining sorted by log-pseudolik; the smaller, the better)
        (A,(C,(B,D)));, with -loglik 30.0
        """)
    nets = readallsnaqnetworks(fn)
    @test length(nets) == 2
    @test SNaQscore.(nets) ≈ [-0.125, -0.3]
    @test SNaQ.loglik.(nets) ≈ [12.5, 30.0]
end
