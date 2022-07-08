module  TestIO

using PopGenCore
using Test

#cats_gen = normpath(joinpath(@__DIR__,"data/", "nancycats.gen"))
#sharks_csv = normpath(joinpath(@__DIR__,"data/", "gulfsharks.csv"))
example_vcf = normpath(joinpath(@__DIR__,"../src/precompile", "precompile.vcf.gz"))
#example_plink = normpath(joinpath(@__DIR__,"../src/precompile", "precompile.vcf.gz"))

@testset "File IO" begin
    @testset "Genepop io" begin
        @test typeof(@nancycats) == PopData
    end

    @testset "delimited io" begin
        @test typeof(@gulfsharks) == PopData
    end

    @testset "VCF io" begin
        @test typeof(vcf(example_vcf, silent = true)) == PopData    
    end

    @testset "Baypass io" begin
        @test_throws ArgumentError baypass(@nancycats)
        sharks = dropmultiallelic(@gulfsharks)
        @test baypass(sharks) isa Matrix
        baypass(sharks, filename = "CI_test.baypass")
        @test isfile("CI_test.baypass")
        rm("CI_test.baypass")
    end
end

end # module

