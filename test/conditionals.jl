module  TestConditionals

using PopGenCore
using Test

x = @nancycats
a = (Int8(1),Int8(2))
b = (Int8(1),Int8(1))
c = (Int8(1),Int8(2), Int8(1))
d = (Int8(2),Int8(2), Int8(2))

@testset "Conditionals.jl" begin
    @testset "biallelic" begin
        @test isbiallelic(x) == false
        @test isbiallelic(fill((1,2), 10)) == true
        @test isbiallelic(x.genodata) == false
    end

    @testset "homozygous" begin
        @test ishom(missing) == false
        @test PopGenCore._ishom(missing) === missing
        @test ishom(a) == false
        @test ishom(b) == true
        @test ishom(c) == false
        @test ishom(d) == true
        @test all(ishom(x.genodata.genotype[1:5]) .== [false, false, false, false, false])
        @test all(ishom(skipmissing(x.genodata.genotype[1:5])) .== [false, false, false])
        @test ishom(b, 1) == true
        @test ishom(a, 1) == false
        @test ishom(b, 2) == false
        @test ishom(missing, 2) == false
        @test PopGenCore._ishom(missing, 2) === missing
    end

    @testset "heterozygous" begin
        @test ishet(missing) == false
        @test PopGenCore._ishet(missing) === missing
        @test ishet(a) == true
        @test ishet(b) == false
        @test ishet(c) == true
        @test ishet(d) == false
        @test all(ishet(x.genodata.genotype[1:5]) .=== [false, false, true, true, true])
        @test all(ishet(skipmissing(x.genodata.genotype[1:5])) .== [true, true, true])
        @test ishet(b, 1) == false
        @test ishet(a, 1) == true
        @test ishet(b, 2) == false
        @test ishet(missing, 2) == false
        @test PopGenCore._ishet(missing, 2) === missing
    end
end

end # module
