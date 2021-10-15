module  TestAlleleGenoFreq

using PopGenCore
using DataFrames
using Test

cats = @nancycats ;
g = cats.genodata.genotype[1:10]
sharks = @gulfsharks ;

@testset "Allele Frequncies" begin
    @testset "Basic" begin
        @test allele_freq(cats) isa NamedTuple
        @test allele_freq(g) isa Dict
        @test allele_freq(g[5]) isa Dict
        @test allele_freq_vec(g) isa AbstractVector
        a = allele_freq(g)
        b = allele_freq(cats.genodata.genotype[11:20])
        @test avg_allele_freq([a,b]) isa Dict
        @test avg_allele_freq([a,b], 2) isa Dict
        @test allele_freq(cats, "fca8") isa Dict
        @test allele_freq(cats, "fca8", population = true) isa DataFrame
        @test allele_freq(133, g) isa AbstractFloat
    end

    @testset "Nuanced" begin
        @test length(allele_freq(cats)) == 9
        @test eltype(allele_freq(cats)) == Dict{Int16,Float64}
        @test allele_freq(cats.genodata.genotype) isa Dict{Int16,Float64}
        df_cats = DataFrames.combine(
            groupby(cats.genodata, :population),
            :genotype => allele_freq => :genos
        )
        @test size(df_cats) == (17,2)
        @test eltype(df_cats.genos) == Dict{Int16,Float64}
        @test length(allele_freq(cats, "fca8")) == 16
        @test allele_freq(cats, "fca8") isa Dict{Int16,Float64}
        
        @test length(allele_freq(sharks)) == 2209
        @test eltype(allele_freq(sharks)) == Dict{Int8,Float64}
        @test allele_freq(sharks.genodata.genotype) isa Dict{Int8,Float64}
        df_sharks = DataFrames.combine(
            groupby(sharks.genodata, :population),
            :genotype => allele_freq => :genos
        )
        @test size(df_sharks) == (7,2)
        @test eltype(df_sharks.genos) == Dict{Int8,Float64}
        @test length(allele_freq(sharks, "contig_2784")) == 2
        @test allele_freq(sharks, "contig_2784") isa Dict{Int8,Float64}

        @test allele_freq_vec([(1,1), (2,2), (2,1)]) == [0.5, 0.5]
        @test allele_freq_vec(missing) === missing
    end
end

@testset "Genotype Frequncies" begin
    @testset "Basic" begin
        @test geno_freq(g) isa Dict
        @test geno_freq(cats, "fca8") isa Dict
        @test geno_freq(cats, "fca8", population = true) isa DataFrame
        @test geno_freq_expected(g) isa Dict
        @test geno_freq_expected(cats, "fca8") isa Dict
        @test geno_freq_expected(cats, "fca8", population = true) isa DataFrame
    end

    @testset "Nuanced" begin
        @test length(geno_freq(cats.genodata.genotype)) == 295
        @test typeof(geno_freq(cats.genodata.genotype)) <: Dict{<:Tuple,Float64}
        @test length(geno_freq(cats, "fca8")) == 51
        @test typeof(geno_freq(cats, "fca8")) <: Dict{<:Tuple,Float64}
        @test size(geno_freq(cats, "fca8", population = true)) == (17,2)

        @test length(geno_freq_expected(cats.genodata.genotype)) == 6241
        @test typeof(geno_freq_expected(cats.genodata.genotype)) <: Dict{<:Tuple,Float64}
        @test length(geno_freq_expected(cats, "fca8")) == 256
        @test typeof(geno_freq_expected(cats, "fca8")) <: Dict{<:Tuple,Float64}
        @test size(geno_freq_expected(cats, "fca8", population = true)) == (17,2)

        @test length(geno_freq(sharks.genodata.genotype)) == 12
        @test length(geno_freq(sharks, "contig_2784")) == 2
        @test typeof(geno_freq(sharks, "contig_2784")) <: Dict{<:Tuple,Float64}
        @test size(geno_freq(sharks, "contig_2784", population = true)) == (7,2)

        @test length(geno_freq_expected(sharks.genodata.genotype)) == 25
        @test typeof(geno_freq_expected(sharks.genodata.genotype)) <: Dict{<:Tuple,Float64}
        @test length(geno_freq_expected(sharks, "contig_2784")) == 4
        @test typeof(geno_freq_expected(sharks, "contig_2784")) <: Dict{<:Tuple,Float64}
        @test size(geno_freq_expected(sharks, "contig_2784", population = true)) == (7,2)
    end
end

@testset "Genotype Counts" begin
    @testset "Basic" begin
        @test geno_count_observed(g) isa Dict
        @test geno_count_observed([missing, missing]) === missing
        @test geno_count_expected(g) isa Dict
    end

    @testset "Nuanced" begin
        @test length(geno_count_observed(cats.genodata.genotype)) == 295
        @test typeof(geno_count_observed(cats.genodata.genotype)) <: Dict{<:Tuple,Int64}
        @test length(geno_count_expected(cats.genodata.genotype)) == 6241
        @test typeof(geno_count_expected(cats.genodata.genotype)) <: Dict{<:Tuple,Float64}

        @test length(geno_count_observed(sharks.genodata.genotype)) == 12
        @test typeof(geno_count_observed(sharks.genodata.genotype)) <: Dict{<:Tuple,Int64}
        @test length(geno_count_expected(sharks.genodata.genotype)) == 25
        @test typeof(geno_count_expected(sharks.genodata.genotype)) <: Dict{<:Tuple,Float64}
    end
end

end