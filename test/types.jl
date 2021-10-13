module  TestTypesDatasets

using PopGenCore
using DataFrames
using PooledArrays
using Test

cats = @nancycats;
sharks = @gulfsharks;

@testset "Data Types and Datasets" begin
    @testset "Dataset dimensions" begin
        @test size(cats.sampleinfo) == (237,3)
        @test size(cats.genodata) == (2133,4)
        @test size(sharks.sampleinfo) == (212,5)
        @test size(sharks.genodata) == (468308,4)
    end

    @testset "Dataset column names" begin
        meta_names_sorted = ["latitude","longitude","name","ploidy","population"]
        loci_names_sorted = ["genotype","locus","name","population"]
        @test sort(names(cats.sampleinfo)) == meta_names_sorted[3:end]
        @test sort(names(cats.genodata)) == loci_names_sorted
        @test sort(names(sharks.sampleinfo)) == meta_names_sorted
        @test sort(names(sharks.genodata)) == loci_names_sorted
    end

    @testset "Nancycats column types" begin
        @test eltype(cats.sampleinfo.name) <: AbstractString
        @test eltype(cats.sampleinfo.population) <: AbstractString
        @test typeof(cats.sampleinfo.ploidy) == Vector{Int8}
        @test_throws ArgumentError typeof(cats.sampleinfo.latitude)
        @test_throws ArgumentError typeof(cats.sampleinfo.longitude)
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
        @test eltype(sharks.sampleinfo.name) <: AbstractString
        @test eltype(sharks.sampleinfo.population) <: AbstractString
        @test eltype(sharks.sampleinfo.ploidy) == Int8
        @test eltype(sharks.sampleinfo.latitude) ==  Float64
        @test eltype(sharks.sampleinfo.longitude) == Float64
        @test typeof(sharks.genodata.name) <: PooledArray
        @test eltype(sharks.genodata.name) <: AbstractString
        @test typeof(sharks.genodata.population) <: PooledArray
        @test eltype(sharks.genodata.population) <: AbstractString
        @test typeof(sharks.genodata.locus) <: PooledArray
        @test eltype(sharks.genodata.locus) <: AbstractString
        @test typeof(sharks.genodata.genotype) <: GenoArray
        @test typeof(sharks.genodata.genotype) <: PooledArray
        @test eltype(sharks.genodata.genotype) <: Union{Missing, Genotype}
    end
end

end  # module
