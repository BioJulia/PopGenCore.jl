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
end

end # module

