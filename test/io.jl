module  TestIO

using PopGenIO
using VariantCallFormat, GZip
using Test

#cats_gen = normpath(joinpath(@__DIR__,"data/", "nancycats.gen"))
#sharks_csv = normpath(joinpath(@__DIR__,"data/", "gulfsharks.csv"))
#example_vcf = normpath(joinpath(@__DIR__,"data/", "example.vcf"))

@testset "Genepop io" begin
    @test typeof(@nancycats) == PopData
end

@testset "delimited io" begin
    @test typeof(@gulfsharks) == PopData
end

@testset "VCF io" begin
    @test typeof(vcf("data/example_vcf", silent = true)) == PopData    
end

end # module

