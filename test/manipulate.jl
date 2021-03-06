module  TestManipulate

using PopGenCore
using DataFrames
using StatsBase
using Test

cats = @nancycats;
sharks = @gulfsharks;

@testset "Manipulate.jl" begin
    @testset "locations" begin
        x = rand(length(samplenames(cats))) ; y = rand(length(samplenames(cats)))
        locationdata!(cats, longitude = x, latitude = y)    
        @test cats.sampleinfo.longitude == x
        @test cats.sampleinfo.latitude == y
        @test locationdata(cats).longitude == x
        @test locationdata(cats).latitude == y
    end

    @testset "decimal-minutes locations" begin
        x = fill("11 22.33W", length(samplenames(cats))) ; y = fill("-41 31.52", length(samplenames(cats)))
        locationdata!(cats, longitude = x, latitude = y)
        @test all(cats.sampleinfo.longitude .== -11.3722)
        @test all(cats.sampleinfo.latitude .== -41.5253)
    end

    @testset "loci and genotypes" begin
        @test length(loci(cats)) == 9
        @test genotypes(cats, "N115" => "fca8") == (135,135)
        N115 = genotypes(cats, "N115")
        @test length(N115) == 9
        @test N115 isa GenoArray
        @test size(genotypes(cats, ["N115", "N7"] => "fca8")) == (2,4)
        @test size(genotypes(cats, "N115" => ["fca8", "fca37"])) == (2,4)
        @test size(genotypes(cats, ["N115", "N7"] => ["fca8", "fca37"])) == (4,4)
        @test length(genotypes(sharks, "contig_475")) == 212
    end

    @testset "populations" begin
        @test length(populations(cats)) == 17
        @test populations(cats, counts = true).population isa Vector{String}
        @test populations(cats, counts = true).count isa Vector{Int}

        rn_dict = Dict("1" => "one", "2" => "two")
        populations!(cats, rn_dict)
        @test "one" ∈ cats.sampleinfo.population && "1" ∉ cats.sampleinfo.population
        @test "two" ∈ cats.sampleinfo.population && "2" ∉ cats.sampleinfo.population

        rn_vect = string.(1:17)
        populations!(cats, rn_vect)
        @test "one" ∉ cats.sampleinfo.population && "1" ∈ cats.sampleinfo.population
        @test "two" ∉ cats.sampleinfo.population && "2" ∈ cats.sampleinfo.population
        rn_vect = string.(1:237)
        populations!(cats, rn_vect)
        @test cats.metadata.populations == 237
        @test cats.sampleinfo.population == rn_vect
        @test_throws DimensionMismatch populations!(cats, string.(1:22)) 

        rn_vect_ii = ["N215", "N297"]
        rn_vect_iin = ["one", "seventeen"]
        populations!(cats, rn_vect_ii, rn_vect_iin)
        @test "one" ∈ cats.sampleinfo.population && "seventeen" ∈ cats.sampleinfo.population
        @test cats.sampleinfo[cats.sampleinfo.name .== "N215", :population] == ["one"]
        @test cats.sampleinfo[cats.sampleinfo.name .== "N297", :population] == ["seventeen"]
    end

    @testset "exclusion" begin
        cats = @nancycats ;
        tmp = exclude(cats, name = "N100", population = ["1", "15"])
        @test length(samplenames(tmp)) == 215
        @test length(populations(tmp)) == 15

        tmp = exclude(cats, name = "N7", locus= "fca8", population = "3")
        @test length(loci(tmp)) == 8
        @test length(populations(tmp)) == 16
        @test length(samplenames(tmp)) == 224

        tmp = keep(cats, name = "N104")
        @test length(samplenames(tmp)) == 1
        tmp = keep(cats, name = ["N104", "N105"])
        @test length(samplenames(tmp)) == 2
        tmp = keep(cats, locus = "fca8")
        @test length(loci(tmp)) == 1
        tmp = keep(cats, locus = ["fca8","fca37"])
        @test length(loci(tmp)) == 2
        tmp = keep(cats, population = "1")
        @test length(populations(tmp)) == 1
        tmp = keep(cats, population = 1:2)
        @test length(populations(tmp)) == 2
    end
end

end # module