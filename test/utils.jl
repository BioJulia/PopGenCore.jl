module  TestUtils

using PopGenCore
using Test

x = @nancycats

@testset "Utilities" begin
    @testset "GeneralUtils.jl" begin
        @test typeof(copy(x)) == PopData
        @test size(x) == (;samples = 237, loci = 9)
        @test sort((2,1,3)) == (1,2,3)
        @test convert_coord("-41 31.52") == -41.5253f0
        @test convert_coord.(["-41:31.52", "25 11:54S"]) == Float32[-41.5253, -25.1983]
        @test length(loci(drop_monomorphic(x))) == 9
        @test length(loci(drop_multiallelic(x))) == 0
        tmp = generate_meta(x.genodata)
        @test names(tmp) == ["name", "population", "ploidy", "longitude", "latitude"]
        @test tmp.name == x.genodata.name.pool
        @test all(.!ismissing(tmp.ploidy))
        @test intersect(x.genodata.population.pool, unique(tmp.population)) == x.genodata.population.pool
        @test loci(x) == x.genodata.locus.pool
        @test samples(x) == x.metadata.name
    end


    @testset "GenotypeUtils.jl" begin
        @test allele_count([(1,2), (1,4)]) == 3
        @test allele_count(x.genodata.genotype[1:10]) == 4
        @test alleles([(1,2), (1,4)]) == [1,2,1,4]
        @test alleles([(1,2), (1,4)], false) == [1,2,1,4] 
        @test alleles([(1,2), (1,4), missing], false) == [1,2,1,4] 
        @test all(alleles([(1,2), (1,4), missing], true) .=== [1,2,1,4,missing])
        @test unique_alleles([(1,2), (1,4), missing]) == [1,2,4]
        @test size(loci_dataframe(x)) == (9,237)
        @test size(loci_matrix(x)) == (237,9)
        @test length(phased_matrix(x)) == 2
        @test size.(phased_matrix(x)) == [(237,9),(237,9)]
        @test eltype(phased_matrix(x)) == Matrix{Union{Missing, Int16}}
        @test eltype(phased_matrix(x)[1]) == Union{Missing, Int16}
    end


    @testset "ioUtils.jl" begin
        @test find_ploidy(x.genodata.genotype[1:20]) == 2
        @test find_ploidy([(1,2,3), (1,2,2)]) == 3
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
        @test count_nonzeros([1,2,0,0,3,4,5]) == 5
        @test reciprocal(0) == 0.0
        @test reciprocal(2) == 1/2
        @test reciprocal_sum([1,2,5,0,0,4]) == (1/1 + 1/2 + 1/5 + 1/4)
    end


    @testset "MissingUtils.jl" begin
        @test nonmissing(fill(missing, 5)) == 0
        @test nonmissing([1,2,3,missing, missing]) == 3
        @test nonmissing(x, "fca8") == 217
        @test nonmissings([1,2,missing, missing, 3], [3,missing, missing, 4, 5]) == [1,5]
    end
end
end