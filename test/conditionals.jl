module  TestConditionals

using PopGenCore
using Test

x = @nancycats

@testset "biallelic" begin
    @test isbiallelic(x) == false
    @test isbiallelic(fill((1,2), 10)) == true
    @test isbiallelic(x.genodata) == false
end

@testset "homozygous" begin
    @test ishom(missing) === missing
    @test ishom((1,2)) == false
    @test ishom((1,1)) == true
    @test ishom((1,1,2)) == false
    @test ishom((2,2,2)) == true
    @test all(ishom(x.genodata.genotype[1:5]) .=== [missing, missing, false, false, false])
    @test all(ishom(skipmissing(x.genodata.genotype[1:5])) .== [false, false, false])
    @test ishom((1,1), 1) == true
    @test ishom((1,2), 1) == false
    @test ishom((1,1), 2) == false
    @test ishom(missing, 2) === missing

end

@testset "heterozygous" begin
    @test ishet(missing) === missing
    @test ishet((1,2)) == true
    @test ishet((1,1)) == false
    @test ishet((1,1,2)) == true
    @test ishet((2,2,2)) == false
    @test all(ishet(x.genodata.genotype[1:5]) .=== [missing, missing, true, true, true])
    @test all(ishet(skipmissing(x.genodata.genotype[1:5])) .== [true, true, true])
    @test ishet((1,1), 1) == false
    @test ishet((1,2), 1) == true
    @test ishet((1,1), 2) == false
    @test ishet(missing, 2) === missing
end

end # module
