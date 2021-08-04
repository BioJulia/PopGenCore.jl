module  TestTypesDatasets

using PopGenCore
using DataFrames
using PooledArrays
using Test

cats = @nancycats;
sharks = @gulfsharks;

@testset "Dataset types" begin
    @test typeof(cats) == PopData
    @test typeof(cats.metadata) == DataFrame
    @test typeof(cats.genotypes) == DataFrame
    @test typeof(sharks) == PopData
    @test typeof(sharks.metadata) == DataFrame
    @test typeof(sharks.genotypes) == DataFrame
end

@testset "Dataset dimensions" begin
    @test size(cats.metadata) == (237,5)
    @test size(cats.genotypes) == (2133,4)
    @test size(sharks.metadata) == (212,5)
    @test size(sharks.genotypes) == (468308,4)
end

@testset "Dataset column names" begin
    meta_names_sorted = ["latitude","longitude","name","ploidy","population"]
    loci_names_sorted = ["genotype","locus","name","population"]
    @test sort(names(cats.metadata)) == meta_names_sorted
    @test sort(names(cats.genotypes)) == loci_names_sorted
    @test sort(names(sharks.metadata)) == meta_names_sorted
    @test sort(names(sharks.genotypes)) == loci_names_sorted
end

@testset "Nancycats column types" begin
    @test typeof(cats.metadata.name) == Vector{String}
    @test typeof(cats.metadata.population) == Vector{String}
    @test typeof(cats.metadata.ploidy) == Vector{Int8}
    @test typeof(cats.metadata.latitude) ==  Vector{Union{Missing, Float32}}
    @test typeof(cats.metadata.longitude) == Vector{Union{Missing, Float32}}
    @test typeof(cats.genotypes.name) <: PooledArray
    @test eltype(cats.genotypes.name) == String
    @test typeof(cats.genotypes.population) <: PooledArray
    @test eltype(cats.genotypes.population) == String
    @test typeof(cats.genotypes.locus) <: PooledArray
    @test eltype(cats.genotypes.locus) == String
    @test typeof(cats.genotypes.genotype) <: GenoArray
    @test eltype(cats.genotypes.genotype) <: Union{Missing, Genotype}
end

@testset "Gulfsharks column types" begin
    @test typeof(sharks.metadata.name) == Vector{String}
    @test typeof(sharks.metadata.population) == Vector{String}
    @test typeof(sharks.metadata.ploidy) == Vector{Int8}
    @test typeof(sharks.metadata.latitude) ==  Vector{Float64}
    @test typeof(sharks.metadata.longitude) == Vector{Float64}
    @test typeof(sharks.genotypes.name) <: PooledArray
    @test eltype(sharks.genotypes.name) == String
    @test typeof(sharks.genotypes.population) <: PooledArray
    @test eltype(sharks.genotypes.population) == String
    @test typeof(sharks.genotypes.locus) <: PooledArray
    @test eltype(sharks.genotypes.locus) == String
    @test typeof(sharks.genotypes.genotype) <: GenoArray
    @test eltype(sharks.genotypes.genotype) <: Union{Missing, Genotype}
end


end  # module
