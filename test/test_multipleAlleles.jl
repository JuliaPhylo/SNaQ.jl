@testset "multiple alleles" begin
import SNaQ: sorttaxa!
global tree, df, d, net, currT

@testset "test: map alleles to species" begin
    tree = readnewick("(6,(5,(7,(3,4))));");
    SNaQ.expandLeaves!(["7"],tree)
    @test writenewick(tree) == "(6,(5,((7:0.0,7__2:0.0):1.0,(3,4))));"
    SNaQ.mergeLeaves!(tree)
    @test writenewick(tree) == "(6,(5,(7:1.0,(3,4))));"
    alleleDF=DataFrame(allele=["1","2"], species=["7","7"])
    CSV.write("tmp.csv", alleleDF);
    df = (@test_logs (:warn, r"^not all alleles were mapped") mapallelesCFtable("tmp.csv",
      joinpath(@__DIR__, "..", "examples", "tableCFCI.csv"),
      # joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "tableCFCI.csv"),
      filename="CFmapped.csv"))
    rm("CFmapped.csv")
    rm("tmp.csv")
    @test df[!,:t4] == ["4","7","3","7","3","3","7","3","3","3","7","3","3","3","3"]
end

#----------------------------------------------------------#
#   testing sorting of taxa and CFs                        #
#----------------------------------------------------------#
@testset "sorttaxa!" begin

letters = ["a","b","c","d"]; cfvalues = [0.6, 0.39, 0.01] # for ab_cd, ac_bd, ad_bc
d = DataFrame(t1=Array{String}(undef,24),t2=Array{String}(undef,24),t3=Array{String}(undef,24),t4=Array{String}(undef,24),
              CF12_34=Array{Float64}(undef,24), CF13_24=Array{Float64}(undef,24), CF14_23=Array{Float64}(undef,24));
irow=1        # d will contain 6!=24 rows: for all permutations on 4 letters
for i1 in 1:4
  ind234 = deleteat!(collect(1:4),i1)
  for i2 in ind234
    ind34 = deepcopy(ind234)
    deleteat!(ind34, findfirst(isequal(i2), ind34))
    for j in 1:2
      i3=ind34[j]; i4=ind34[3-j]
      d[irow,:t1]=letters[i1]; d[irow,:t2]=letters[i2]; d[irow,:t3]=letters[i3]; d[irow,:t4]=letters[i4]
      # CF12_34 corresponds to CFi1i2_i3i4
      if     (i1,i2)∈[(1,2),(2,1),(3,4),(4,3)] d[irow,:CF12_34] = cfvalues[1]
      elseif (i1,i2)∈[(1,3),(3,1),(2,4),(4,2)] d[irow,:CF12_34] = cfvalues[2]
      elseif (i1,i2)∈[(1,4),(4,1),(2,3),(3,2)] d[irow,:CF12_34] = cfvalues[3]
      end # next: set CF13_24
      if     (i1,i3)∈[(1,2),(2,1),(3,4),(4,3)] d[irow,:CF13_24] = cfvalues[1]
      elseif (i1,i3)∈[(1,3),(3,1),(2,4),(4,2)] d[irow,:CF13_24] = cfvalues[2]
      elseif (i1,i3)∈[(1,4),(4,1),(2,3),(3,2)] d[irow,:CF13_24] = cfvalues[3]
      end # nest: set CF14_23
      if     (i1,i4)∈[(1,2),(2,1),(3,4),(4,3)] d[irow,:CF14_23] = cfvalues[1]
      elseif (i1,i4)∈[(1,3),(3,1),(2,4),(4,2)] d[irow,:CF14_23] = cfvalues[2]
      elseif (i1,i4)∈[(1,4),(4,1),(2,3),(3,2)] d[irow,:CF14_23] = cfvalues[3]
      end
      irow += 1
    end
  end
end
# d
d2 = deepcopy(d);
sorttaxa!(d2);
d3 = DataFrame(t1=repeat([letters[1]],outer=[24]),t2=repeat([letters[2]],outer=[24]),
               t3=repeat([letters[3]],outer=[24]),t4=repeat([letters[4]],outer=[24]),
               CF12_34=repeat([cfvalues[1]],outer=[24]),CF13_24=repeat([cfvalues[2]],outer=[24]),CF14_23=repeat([cfvalues[3]],outer=[24]));
@test d2==d3

dat = readtableCF(d);
net = (@test_logs readnewick("(a,((b)#H1,((#H1,c),d)));"));
# earlier warning: "net does not have identifiable branch lengths"
sorttaxa!(dat)

@test [q.obsCF for q in dat.quartet] == [[0.6,0.39,0.01] for i in 1:24]
@test [q.taxon for q in dat.quartet] == [letters for i in 1:24]

end # of testset: sorttaxa!

@testset "snaq on multiple alleles" begin

df = DataFrame(t1=["6","7"], t2=["7","6"], t3=["4","4"], t4=["8","8"],
    a=[true,true], # to test recognition of columns
    CF12_34=[0.25, 0.15], ngenes=[10,20],
    CF13_24=[0.3,0.55], b=[false,false], CF14_23=[0.45,0.3])
@test length(readtableCF(df).quartet) == 2
d = readtableCF(df, mergerows=true)
@test isempty(d.repSpecies)
@test length(d.quartet) == 1
@test d.quartet[1].obsCF ≈ [0.3, 0.5, 0.2]
@test d.quartet[1].ngenes ≈ 15
SNaQ.descData(d, devnull)
SNaQ.descData(d, "tmp.log")
summarizedataCF(d, filename="tmp.log")
rm("tmp.log")


# New testing example taken directly from the wiki
mappingfile = joinpath(dirname(pathof(SNaQ)), "..","examples",
           "mappingIndividuals.csv");
tm = CSV.read(mappingfile, DataFrame)
taxonmap = Dict(r[:individual] => r[:species] for r in eachrow(tm))
genetreefile = joinpath(dirname(pathof(SNaQ)), "..","examples",
           "genetrees_alleletips.tre");
genetrees = readmultinewick(genetreefile);
df_sp = tablequartetCF(countquartetsintrees(genetrees, taxonmap;
           showprogressbar=false)...);
d_sp = readtableCF(DataFrame(df_sp))
T_sp = readnewick("((S4,S5),((S1,S3),S2));")
for E in T_sp.edge E.length = 0.25 end
snaq!(T_sp, d_sp)


end # test of snaq on multiple alleles

end # overall multiple allele sets of testests
