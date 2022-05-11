module  TestAlleleGenoFreq

using PopGenCore
using DataFrames
using Test

cats = @nancycats ;
g = cats.genodata.genotype[1:10]
sharks = @gulfsharks ;

@testset "Allele Frequncies" begin
    @testset "Basic" begin
        @test allelefreq(cats) isa NamedTuple
        @test allelefreq(g) isa Dict
        @test allelefreq(g[5]) isa Dict
        @test allelefreq_vec(g) isa AbstractVector
        a = allelefreq(g)
        b = allelefreq(cats.genodata.genotype[11:20])
        @test avg_allelefreq([a,b]) isa Dict
        @test avg_allelefreq([a,b], 2) isa Dict
        @test allelefreq(cats, "fca8") isa Dict
        @test allelefreq(cats, "fca8", population = true) isa DataFrame
        @test allelefreq(g, 133) isa AbstractFloat
    end

    @testset "Nuanced" begin
        @test length(allelefreq(cats)) == 9
        @test eltype(allelefreq(cats)) == Dict{Int16,Float64}
        @test allelefreq(cats.genodata.genotype) isa Dict{Int16,Float64}
        df_cats = DataFrames.combine(
            groupby(cats.genodata, :population),
            :genotype => allelefreq => :genos
        )
        @test size(df_cats) == (17,2)
        @test eltype(df_cats.genos) == Dict{Int16,Float64}
        @test length(allelefreq(cats, "fca8")) == 16
        @test allelefreq(cats, "fca8") isa Dict{Int16,Float64}
        
        @test length(allelefreq(sharks)) == 2209
        @test eltype(allelefreq(sharks)) == Dict{Int8,Float64}
        @test allelefreq(sharks.genodata.genotype) isa Dict{Int8,Float64}
        df_sharks = DataFrames.combine(
            groupby(sharks.genodata, :population),
            :genotype => allelefreq => :genos
        )
        @test size(df_sharks) == (7,2)
        @test eltype(df_sharks.genos) == Dict{Int8,Float64}
        @test length(allelefreq(sharks, "contig_2784")) == 2
        @test allelefreq(sharks, "contig_2784") isa Dict{Int8,Float64}

        @test allelefreq_vec([(1,1), (2,2), (2,1)]) == [0.5, 0.5]
        @test allelefreq_vec(missing) === missing
    end
end

@testset "Genotype Frequncies" begin
    @testset "Basic" begin
        @test genofreq(g) isa Dict
        @test genofreq(cats, "fca8") isa Dict
        @test genofreq(cats, "fca8", population = true) isa DataFrame
        @test genofreq_expected(g) isa Dict
        @test genofreq_expected(cats, "fca8") isa Dict
        @test genofreq_expected(cats, "fca8", population = true) isa DataFrame
    end

    @testset "Nuanced" begin
        @test length(genofreq(cats.genodata.genotype)) == 295
        @test typeof(genofreq(cats.genodata.genotype)) <: Dict{<:Tuple,Float64}
        @test length(genofreq(cats, "fca8")) == 51
        @test typeof(genofreq(cats, "fca8")) <: Dict{<:Tuple,Float64}
        @test size(genofreq(cats, "fca8", population = true)) == (17,2)

        @test length(genofreq_expected(cats.genodata.genotype)) == 3160
        @test typeof(genofreq_expected(cats.genodata.genotype)) <: Dict{<:Tuple,Float64}
        @test length(genofreq_expected(cats, "fca8")) == 136
        @test typeof(genofreq_expected(cats, "fca8")) <: Dict{<:Tuple,Float64}
        @test size(genofreq_expected(cats, "fca8", population = true)) == (17,2)

        @test length(genofreq(sharks.genodata.genotype)) == 12
        @test length(genofreq(sharks, "contig_2784")) == 2
        @test typeof(genofreq(sharks, "contig_2784")) <: Dict{<:Tuple,Float64}
        @test size(genofreq(sharks, "contig_2784", population = true)) == (7,2)

        @test length(genofreq_expected(sharks.genodata.genotype)) == 15
        @test typeof(genofreq_expected(sharks.genodata.genotype)) <: Dict{<:Tuple,Float64}
        @test length(genofreq_expected(sharks, "contig_2784")) == 3
        @test typeof(genofreq_expected(sharks, "contig_2784")) <: Dict{<:Tuple,Float64}
        @test size(genofreq_expected(sharks, "contig_2784", population = true)) == (7,2)
    end
end

@testset "Genotype Counts" begin
    @testset "Basic" begin
        @test genocount_observed(g) isa Dict
        @test genocount_observed([missing, missing]) === missing
        @test genocount_expected(g) isa Dict
    end

    @testset "Nuanced" begin
        @test length(genocount_observed(cats.genodata.genotype)) == 295
        @test typeof(genocount_observed(cats.genodata.genotype)) <: Dict{<:Tuple,Int64}
        @test length(genocount_expected(cats.genodata.genotype)) == 3160
        @test typeof(genocount_expected(cats.genodata.genotype)) <: Dict{<:Tuple,Float64}

        @test length(genocount_observed(sharks.genodata.genotype)) == 12
        @test typeof(genocount_observed(sharks.genodata.genotype)) <: Dict{<:Tuple,Int64}
        @test length(genocount_expected(sharks.genodata.genotype)) == 15
        @test typeof(genocount_expected(sharks.genodata.genotype)) <: Dict{<:Tuple,Float64}
    end
end

end