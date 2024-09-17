using SNaQ
using Test
using Aqua

@testset "SNaQ.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SNaQ)
    end
    # Write your tests here.
end
