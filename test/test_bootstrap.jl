# test the functions in src/bootstrap.jl

exdir = joinpath(@__DIR__,"..","examples")
# exdir = joinpath(dirname(pathof(PhyloNetworks)), "..","examples")



@testset "bootsnaq from quartet CF intervals" begin
T=readnewick(joinpath(exdir,"startTree.txt"))
datf=DataFrame(CSV.File(joinpath(exdir,"tableCFCI.csv")); copycols=false)
originalstdout = stdout
redirect_stdout(devnull) # requires julia v1.6
bootnet = bootsnaq(T,datf,nrep=2,runs=1,seed=1234,filename="",Nfail=2,
                   ftolAbs=1e-3,ftolRel=1e-3,xtolAbs=1e-4,xtolRel=1e-3,liktolAbs=0.01)
redirect_stdout(originalstdout)
@test size(bootnet)==(2,)
@test writeTopology(bootnet[1], round=true) != writeTopology(bootnet[2], round=true)
# "(2,((5,#H9:0.0::0.298):3.927,3):1.331,(((1,6):0.019,4):0.0)#H9:0.0::0.702);"
# above: bad diamond 2, and both edges above the hybrid have estimated length of 0.0...
end


@testset "bootsnaq from bootstrap gene trees, multiple procs" begin
treefile = joinpath(exdir,"treefile.txt") # pretending these are bootstrap trees, for all genes
boottrees = Vector{HybridNetwork}[]
T=readnewick(joinpath(exdir,"startTree.txt"))
net1 = readnewick("(5,(((2,(1)#H7:::0.7143969494428192):1.5121337017411736,4):0.4894187322508883,3):0.519160762355313,(6,#H7:::0.2856030505571808):10.0);")
for i=1:13 push!(boottrees, readmultinewick(treefile)) end
for i=1:13 @test size(boottrees[i])==(10,) end # 10 bootstrap trees for each of 13 "genes"
addprocs(1)
@everywhere using SNaQ
# using Distributed; @everywhere begin; using Pkg; Pkg.activate("."); using PhyloNetworks; end
originalstdout = stdout
redirect_stdout(devnull)
bootnet = bootsnaq(T,boottrees,nrep=2,runs=2,otherNet=net1,seed=1234,
                   prcnet=0.5,filename="",Nfail=2,ftolAbs=1e-3,ftolRel=1e-3)
redirect_stdout(originalstdout)
rmprocs(workers())
@test size(bootnet)==(2,)
@test all(n -> n.numHybrids==1, bootnet)
@test writeTopology(bootnet[1], round=true, digits=1) != writeTopology(bootnet[2], round=true, digits=1)
filelist = joinpath(exdir, "treefilelist.txt")
boottrees = readBootstrapTrees(filelist)
@test length(boottrees) == 2
@test [length(b) for b in boottrees] == [10, 10]
end
