## Claudia: commented these tests out because they were not part of runtests.jl
## when I added the test for readNexusTrees. So, I was not sure if they were
## not included for a reason (slow test?)

@testset "Read CF data" begin
  examp_dir = joinpath(@__DIR__, "..", "examples")
  d = readtableCF(joinpath(examp_dir, "tableCF.txt"))
  d = readtableCF(joinpath(examp_dir, "tableCF_species.csv"))
  d = readtableCF(joinpath(examp_dir, "tableCFCI.csv"))

  sixtreestr = ["(E,((A,B),(C,D)),O);","(((A,B),(C,D)),(E,O));","(A,B,((C,D),(E,O)));",
              "(B,((C,D),(E,O)));","((C,D),(A,(B,E)),O);","((C,D),(A,B,E),O);"]
  sixtrees = readnewick.(sixtreestr)
  d = DataFrame(tablequartetCF(readtrees2CF(sixtrees, writeTab=false, writeSummary=false)))
  @test nrow(d) > 0
  d = DataFrame(tablequartetCF(readtrees2CF(sixtrees, writeTab=false, writeSummary=false, whichQ="rand")))
  @test nrow(d) > 0

  for nq in [1, 3, 7, 12]
    d = DataFrame(tablequartetCF(readtrees2CF(sixtrees, writeTab=false, writeSummary=false, whichQ="rand", numQ=nq)))
    @test nrow(d) == nq
  end

  d = readtrees2CF(sixtrees)
  SNaQ.descData(d)
end


#Test consistency of writing/reading. Moved from test_relaxed_reading. Now only in PN. 
net = readnewicklevel1("(E,((B)#H1:::.5,((D,C),(#H1:::.5,A))));");
@test writenewick(net) == "(D:1.0,C:1.0,((#H1:1.0::0.5,A:1.0):1.0,((B:1.0)#H1:1.0::0.5,E:1.0):1.0):1.0);"
originalstdout = stdout
redirect_stdout(devnull) # requires julia v1.6
@test_logs SNaQ.printEverything(net)
redirect_stdout(originalstdout)

@testset "test: reading nexus file" begin
# with translate table, hybrids, failed γ in net 2, bad net 3, v2 format in net 4
nexusfile = joinpath(@__DIR__, "..", "examples", "test_reticulatetreeblock.nex")
# nexusfile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","test_reticulatetreeblock.nex")
vnet = (@test_logs (:warn, r"^hybrid edge") (:warn,r"^skipped") readnexus_treeblock(nexusfile));
@test length(vnet) == 3
@test writenewick(vnet[1]) == "((tax4,(tax3,#H7:0.001::0.08):0.3):0.6,(tax2,(tax1:0.1)#H7:0.9::0.92):10.0);"
@test writenewick(vnet[3]) == "((tax1:1.13,((tax2:0.21)#H1:0.89::0.72,(tax3:1.03,(#H1:0.3::0.28,tax4:0.51)S3:0.51)S4:0.08)S5:0.2):0.6,tax5:1.14);"
# example without translate table and without reticulations
nexusfile = joinpath(@__DIR__, "..", "examples", "test.nex")
# nexusfile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","test.nex")
vnet = readnexus_treeblock(nexusfile, SNaQ.readTopologyUpdate, false, false; reticulate=false);
@test length(vnet) == 10
@test length(vnet[10].edge) == 9
@test vnet[10].edge[7].length ≈ 0.00035
end


@testset "test: calculate quartet CF from input gene trees" begin
sixtreestr = ["(E,((A,B),(C,D)),O);","(((A,B),(C,D)),(E,O));","(A,B,((C,D),(E,O)));",
              "(B,((C,D),(E,O)));","((C,D),(A,(B,E)),O);","((C,D),(A,B,E),O);"]
sixtrees = readnewick.(sixtreestr)
df1 = DataFrame(tablequartetCF(countquartetsintrees(sixtrees)...))
df2 = DataFrame(tablequartetCF(readtrees2CF(sixtrees, writeTab=false, writeSummary=false)))
o = [1,2,4,7,11,3,5,8,12,6,9,13,10,14,15]
select!(df1, Not(:qind))
@test df1 == df2[o,:]
## weight each allele, countquartetsintrees: already tested in PN
# q,t = countquartetsintrees(sixtrees, Dict("A"=>"AB", "B"=>"AB"); weight_byallele=true, showprogressbar=false);
# df1 = tablequartetCF(q,t)
# mapallelesCFtable performs a different averaging: first across each set of 4 alleles,
# then across sets of 4 alleles that map to the same set of 4 species
CSV.write("tmp_qCF.csv", df2)
CSV.write("tmp_map.csv", DataFrame(allele = ["A","B"], species = ["AB","AB"]))
df2_byallele = (@test_logs (:warn, r"not all alleles were mapped") mapallelesCFtable("tmp_map.csv", "tmp_qCF.csv"))
rm.(["tmp_qCF.csv","tmp_map.csv"]);
q = readtableCF!(df2_byallele);
df2 = tablequartetCF(q) # 45×8 NamedTuple
# df2[7:11,:]
# df12 = innerjoin(df1, df2, on=[:t1,:t2,:t3,:t4], makeunique=true)
# all([df12[:,4+i] ≈ df12[:,8+i] for i in 1:4]) # false: because averaging done differently by the 2 functions
@test df2 == (
  # add qind?
  t1=["AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","C"],
  t2=["AB__2","AB__2","AB__2","AB__2","AB__2","AB__2","C","C","C","D","D"],
  t3=["C","C","C","D","D","E","D","D","E","E","E"],
  t4=["D","E","O","E","O","O","E","O","O","O","O"],
  CF12_34=[1.0,0.75,1.0,0.75,1.0,0.75,0.0,0.0,(3/5+4/6)/2,(3/5+4/6)/2,1.0],
  CF13_24=[0,0.25,0,0.25,0,0,0,0,(2/5+2/6)/2,(2/5+2/6)/2,0.0],
  CF14_23=[0,0,0,0,0,0.25,1,1,0,0,0],
  ngenes=Union{Missing,Float64}[5,4,5,4,5,4,5.5,5.5,5.5,5.5,6]
)
end