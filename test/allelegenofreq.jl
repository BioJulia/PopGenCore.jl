module  TestAlleleGenoFreq

using PopGenCore
using Test

cats = @nancycats
g = cats.genodata.genotype[1:10]

@testset "Allele Frequencies" begin
    @test allele_freq(cats) isa NamedTuple
    @test allele_freq(g) isa Dict
    @test allele_freq(g[5]) isa Dict
    @test allele_freq_vec(g) isa AbstractVector
    a = allele_freq(g)
    b = allele_freq(cats.genodata.genotype[11:20])
    @test avg_allele_freq([a,b]) isa Dict
    @test avg_allele_freq([a,b], 2) isa Dict
    @test allele_freq(cats, "fca8") isa Dict
    @test allele_freq(cats, "fca8", population = true) isa DataFrames.DataFrame
    @test allele_freq(133, g) isa AbstractFloat
end

@testset "Genotype Frequencies" begin
    @test geno_count_observed(g) isa Dict
    @test geno_count_observed([missing, missing]) === missing
    @test geno_count_expected(g) isa Dict
    @test geno_freq(g) isa Dict
    @test geno_freq(cats, "fca8") isa Dict
    @test geno_freq(cats, "fca8", population = true) isa DataFrames.DataFrame
    @test geno_freq_expected(g) isa Dict
    @test geno_freq_expected(cats, "fca8") isa Dict
    @test geno_freq_expected(cats, "fca8", population = true) isa DataFrames.DataFrame
end

end