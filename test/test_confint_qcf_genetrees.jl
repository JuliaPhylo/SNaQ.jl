# let us make a toy example with the file examples/raxmltrees.tre
# we will have thre gene trees and three their bootstrap samples, each with 10 replicates
# raxmltrees.tre has 30 trees. Make packages of 10, and pick one at random to be the mle
# then we have an mle and the replicates for a dataset with 3 gene trees

examp_dir = joinpath(@__DIR__, "..", "examples")
trees = readmultinewick(joinpath(examp_dir, "raxmltrees.tre"))

gts_boots = [trees[1:10], trees[11:20], trees[21:30]]

gts_mles = [trees[2], trees[12], trees[22]]

qcf_cis = confintqCF_genetrees(gts_mles, gts_boots, 0.95, false)

@testset "Testing the order of CI and MLE" begin
    @test qcf_cis.CF12_34_lo[1] <= qcf_cis.CF12_34[1] <= qcf_cis.CF12_34_hi[1]
    @test qcf_cis.CF12_34_lo[2] <= qcf_cis.CF12_34[2] <= qcf_cis.CF12_34_hi[2]
    @test qcf_cis.CF12_34_lo[3] <= qcf_cis.CF12_34[3] <= qcf_cis.CF12_34_hi[3]
    @test qcf_cis.CF12_34_lo[4] <= qcf_cis.CF12_34[4] <= qcf_cis.CF12_34_hi[4]
    @test qcf_cis.CF12_34_lo[5] <= qcf_cis.CF12_34[5] <= qcf_cis.CF12_34_hi[5]    
    @test qcf_cis.CF12_34_lo[6] <= qcf_cis.CF12_34[6] <= qcf_cis.CF12_34_hi[6]    
    @test qcf_cis.CF12_34_lo[7] <= qcf_cis.CF12_34[7] <= qcf_cis.CF12_34_hi[7]
    @test qcf_cis.CF12_34_lo[8] <= qcf_cis.CF12_34[8] <= qcf_cis.CF12_34_hi[8]
    @test qcf_cis.CF12_34_lo[9] <= qcf_cis.CF12_34[9] <= qcf_cis.CF12_34_hi[9]
    @test qcf_cis.CF12_34_lo[10] <= qcf_cis.CF12_34[10] <= qcf_cis.CF12_34_hi[10]
    @test qcf_cis.CF12_34_lo[11] <= qcf_cis.CF12_34[11] <= qcf_cis.CF12_34_hi[11]
    @test qcf_cis.CF12_34_lo[12] <= qcf_cis.CF12_34[12] <= qcf_cis.CF12_34_hi[12]
    @test qcf_cis.CF12_34_lo[13] <= qcf_cis.CF12_34[13] <= qcf_cis.CF12_34_hi[13]
    @test qcf_cis.CF12_34_lo[14] <= qcf_cis.CF12_34[14] <= qcf_cis.CF12_34_hi[14]
    @test qcf_cis.CF12_34_lo[15] <= qcf_cis.CF12_34[15] <= qcf_cis.CF12_34_hi[15]
end

@testset "The dimensions of qcf_cis should be 15 quartet rows x 4 taxon cols, 9 qcf cols, and ngene = 14" begin
    nrow, ncol = size(qcf_cis)
    @test nrow == 15
    @test ncol == 14
end

@testset "The type of the output by confint_qCF_genetrees should be a DataFrame or a subtype of it" begin
    @test typeof(qcf_cis) <: DataFrame
end
