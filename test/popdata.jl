module  TestPopData

using PopGenCore
using DataFrames
using Test

x = @nancycats

@testset "PopData" begin
    @testset "Struct types" begin
        @test x isa PopData
        @test x.metadata isa PopDataInfo
        @test x.genodata isa DataFrame
        @test x.sampleinfo isa DataFrame
        @test x.locusinfo isa DataFrame
    end

    @testset "PopDataInfo" begin
        @test x.metadata == x.info
        @test x.sampleinfo == x.metadata.sampleinfo
        @test PopData(x.genodata) isa PopData
        @test x.metadata.locusinfo == x.locusinfo
        @test x.metadata.ploidy == 2
        @test x.metadata.loci == 9
        @test x.metadata.samples == 237
        @test x.metadata.populations == 17
        @test x.metadata.biallelic == false
    end

    @testset "Indexing DataFrame" begin
        @test x[x.genodata.locus .== "fca8", :] isa DataFrame
        @test x[x.genodata.name .== "N100", :] isa DataFrame
        @test x[x.genodata.name .∈ Ref(["N100", "N217"]), :] isa DataFrame
        @test x[x.genodata.locus .∈ Ref(["fca8", "fca37"]), :] isa DataFrame
        @test x[x.genodata.locus .== "fca8", :name] isa AbstractArray
        @test x[(x.genodata.locus .== "fca8") .& (ishom(x.genodata.genotype, 139)), :name] isa AbstractArray
        
    end
    @testset "Indexing PopData" begin
        @test x[x.genodata.locus .== "fca8"] isa PopData
        @test x[x.genodata.name .== "N100"] isa PopData
        @test x[x.genodata.name .∈ Ref(["N100", "N217"])] isa PopData
        @test x[x.genodata.locus .∈ Ref(["fca8", "fca37"])] isa PopData
        @test x[(x.genodata.locus .== "fca8") .& (ishom(x.genodata.genotype, 139))] isa PopData
    end
end

end #module 