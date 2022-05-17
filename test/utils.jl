module  TestUtils

using PopGenCore
using Test

x = @nancycats

@testset "Utilities" begin
    @testset "GeneralUtils.jl" begin
        @test typeof(copy(x)) == PopData
        @test size(x) == (;samples = 237, loci = 9)
        @test sort((2,1,3)) == (1,2,3)
        @test convertcoord("-41 31.52") == -41.5253
        @test convertcoord.(["-41:31.52", "25 11:54S"]) == [-41.5253, -25.1983]
        @test length(loci(dropmonomorphic(x))) == 9
        @test length(loci(dropmultiallelic(x))) == 0
    end


    @testset "GenotypeUtils.jl" begin
        geno1 = Tuple(Int8[1,2])
        geno2 = Tuple(Int8[1,4])
        @test allelecount([(Int8(1),Int8(2)), (Int8(1),Int8(4))]) == 3
        @test allelecount(x.genodata.genotype[1:10]) == 4
        @test alleles([geno1, geno2]) == [1,2,1,4]
        @test alleles([geno1, geno2], false) == [1,2,1,4] 
        @test alleles([geno1, geno2, missing], false) == [1,2,1,4] 
        @test all(alleles([geno1, geno2, missing], true) .=== [1,2,1,4,missing])
        @test uniquealleles([geno1, geno2, missing]) == [1,2,4]
        @test size(locidataframe(x)) == (9,237)
        @test size(locimatrix(x)) == (237,9)
        @test length(phasedmatrix(x)) == 2
        @test size.(phasedmatrix(x)) == [(237,9),(237,9)]
        @test eltype(phasedmatrix(x)) == Matrix{Union{Missing, Int16}}
        @test eltype(phasedmatrix(x)[1]) == Union{Missing, Int16}
    end


    @testset "ioUtils.jl" begin
        @test findploidy(x.genodata.genotype[1:20]) == 2
        @test findploidy([(1,2,3), (1,2,2)]) == 3
        @test phase("128114", Int16, 3) == (114,128)
        @test phase("128114", Int16, 2) == (12,14,81)
        @test phase("128", Int16, 3) == (128,)
        @test eltype(phase("128", Int16, 3)) == Int16
        @test eltype(phase("110", Int8, 3)) == Int8
        @test phase(128114, Int16, 3) == (114,128)
        @test phase(128114, Int16, 2) == (12,14,81)
        @test phase(128, Int16, 3) == (128,)
        @test eltype(phase(128, Int16, 3)) == Int16
        @test eltype(phase(110, Int8, 3)) == Int8
        @test phase(missing, Int8, 124) === missing
        @test unphase((1,2,3,4), digits = 3) == "001002003004"
        @test unphase(missing, digits = 2, ploidy = 2, miss = -9) == "-9"
        @test unphase(missing, digits = 2, ploidy = 2, miss = 0) == "0000"
    end


    @testset "MathUtils.jl" begin
        @test countnonzeros([1,2,0,0,3,4,5]) == 5
        @test reciprocal(0) == 0.0
        @test reciprocal(2) == 1/2
        @test reciprocalsum([1,2,5,0,0,4]) == (1/1 + 1/2 + 1/5 + 1/4)
    end


    @testset "MissingUtils.jl" begin
        @test nonmissing(fill(missing, 5)) == 0
        @test nonmissing([1,2,3,missing, missing]) == 3
        @test nonmissing(x, "fca8") == 217
        @test nonmissings([1,2,missing, missing, 3], [3,missing, missing, 4, 5]) == [1,5]
        @test isallmissing([missing, missing, missing]) == true
        @test isallmissing([1,missing]) == false
        @test isallmissing([1,2]) == false
    end
end
end