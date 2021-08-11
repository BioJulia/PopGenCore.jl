module  TestManipulate

using PopGenCore
using DataFrames
using StatsBase
using Test

cats = @nancycats;
sharks = @gulfsharks;

@testset "locations" begin
    x = rand(length(samples(cats))) ; y = rand(length(samples(cats)))
    locations!(cats, long = x, lat = y)    
    @test cats.metadata.longitude == Float32.(x)
    @test cats.metadata.latitude == Float32.(y)
    @test locations(cats).longitude == Float32.(x)
    @test locations(cats).latitude == Float32.(y)
end

@testset "decimal-minutes locations" begin
    x = fill("11 22.33W", length(samples(cats))) ; y = fill("-41 31.52", length(samples(cats)))
    locations!(cats, long = x, lat = y)
    @test all(cats.metadata.longitude .== Float32(-11.3722))
    @test all(cats.metadata.latitude .== Float32(-41.5253))
end

@testset "loci and genotypes" begin
    @test length(loci(cats)) == 9
    @test get_genotype(cats, sample = "N115", locus = "fca8") == (135,135)
    N115 = get_genotypes(cats, "N115")
    @test length(N115) == 9
    @test typeof(N115) == Vector{Union{Missing, Tuple{Int16,Int16}}}
    @test typeof(get_genotypes(cats, name = "N115" , locus = "fca8")) <: SubDataFrame
    @test names(get_genotypes(cats, name = "N115" , locus = "fca8")) == ["name", "population", "locus", "genotype"]
    @test size(get_genotypes(cats, name = ["N115", "N7"] , locus = "fca8")) == (2,4)
    @test size(get_genotypes(cats, name = "N115" , locus = ["fca8", "fca37"])) == (2,4)
    @test size(get_genotypes(cats, name = ["N115", "N7"] , locus = ["fca8", "fca37"])) == (4,4)
    @test length(genotypes(sharks, "contig_475")) == 212
end

@testset "populations" begin
    @test length(populations(cats)) == 17
    @test typeof(populations(cats, counts = true).population) == Vector{String}
    @test typeof(populations(cats, counts = true).count) == Vector{Int}

    rn_dict = Dict("1" => "one", "2" => "two")
    populations!(cats, rn_dict)
    @test "one" ∈ cats.metadata.population && "1" ∉ cats.metadata.population
    @test "two" ∈ cats.metadata.population && "2" ∉ cats.metadata.population

    rn_vect = string.(1:17)
    populations!(cats, rn_vect)
    @test "one" ∉ cats.metadata.population && "1" ∈ cats.metadata.population
    @test "two" ∉ cats.metadata.population && "2" ∈ cats.metadata.population
 
    rn_vect_ii = ["N215", "N297"]
    rn_vect_iin = ["one", "seventeen"]
    populations!(cats, rn_vect_ii, rn_vect_iin)
    @test "one" ∈ cats.metadata.population && "seventeen" ∈ cats.metadata.population
    @test cats.metadata[cats.metadata.name .== "N215", :population] == ["one"]
    @test cats.metadata[cats.metadata.name .== "N297", :population] == ["seventeen"]
end

@testset "exclusion" begin
    cats = @nancycats ;
    tmp = exclude(cats, name = "N100", population = ["1", "15"])
    @test length(samples(tmp)) == 215
    @test length(populations(tmp)) == 15

    tmp = exclude(cats, name = "N102", locus= "fca8", population = "3")
    @test length(loci(tmp)) == 8
    @test length(populations(tmp)) == 16
    @test length(samples(tmp)) == 225

    tmp = keep(cats, name = "N104")
    @test length(samples(tmp)) == 1
    tmp = keep(cats, name = ["N104", "N105"])
    @test length(samples(tmp)) == 2
    tmp = keep(cats, locus = "fca8")
    @test length(loci(tmp)) == 1
    tmp = keep(cats, locus = "fca37")
    @test length(loci(tmp)) == 2
end

end # module