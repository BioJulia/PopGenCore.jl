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
    @test typeof(cats.genodata) == DataFrame
    @test typeof(sharks) == PopData
    @test typeof(sharks.metadata) == DataFrame
    @test typeof(sharks.genodata) == DataFrame
end

@testset "Dataset dimensions" begin
    @test size(cats.metadata) == (237,5)
    @test size(cats.genodata) == (2133,4)
    @test size(sharks.metadata) == (212,5)
    @test size(sharks.genodata) == (468308,4)
end

@testset "Dataset column names" begin
    meta_names_sorted = ["latitude","longitude","name","ploidy","population"]
    loci_names_sorted = ["genotype","locus","name","population"]
    @test sort(names(cats.metadata)) == meta_names_sorted
    @test sort(names(cats.genodata)) == loci_names_sorted
    @test sort(names(sharks.metadata)) == meta_names_sorted
    @test sort(names(sharks.genodata)) == loci_names_sorted
end

@testset "Nancycats column types" begin
    @test typeof(cats.metadata.name) <: Vector{AbstractString}
    @test typeof(cats.metadata.population) <: Vector{AbstractString}
    @test typeof(cats.metadata.ploidy) == Vector{Int8}
    @test typeof(cats.metadata.latitude) ==  Vector{Union{Missing, Float32}}
    @test typeof(cats.metadata.longitude) == Vector{Union{Missing, Float32}}
    @test typeof(cats.genodata.name) <: PooledArray
    @test eltype(cats.genodata.name) <: AbstractString
    @test typeof(cats.genodata.population) <: PooledArray
    @test eltype(cats.genodata.population) <: AbstractString
    @test typeof(cats.genodata.locus) <: PooledArray
    @test eltype(cats.genodata.locus) <: AbstractString
    @test typeof(cats.genodata.genotype) <: GenoArray
    @test eltype(cats.genodata.genotype) <: Union{Missing, Genotype}
end

@testset "Gulfsharks column types" begin
    @test typeof(sharks.metadata.name) <: Vector{AbstractString}
    @test typeof(sharks.metadata.population) <: Vector{AbstractString}
    @test typeof(sharks.metadata.ploidy) == Vector{Int8}
    @test typeof(sharks.metadata.latitude) ==  Vector{Float64}
    @test typeof(sharks.metadata.longitude) == Vector{Float64}
    @test typeof(sharks.genodata.name) <: PooledArray
    @test eltype(sharks.genodata.name) <: AbstractString
    @test typeof(sharks.genodata.population) <: PooledArray
    @test eltype(sharks.genodata.population) <: AbstractString
    @test typeof(sharks.genodata.locus) <: PooledArray
    @test eltype(sharks.genodata.locus) <: AbstractString
    @test typeof(sharks.genodata.genotype) <: GenoArray
    @test eltype(sharks.genodata.genotype) <: Union{Missing, Genotype}
end


end  # module
