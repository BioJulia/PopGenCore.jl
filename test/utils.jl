module  TestUtils

using PopGenCore
using Test

x = @nancycats
#=
@testset "general utils" begin

end

@testset "genotype utils" begin

end


@testset "io utils" begin


end
=#
@testset "math utils" begin
    @test count_nonzeros([1,2,0,0,3,4,5]) == 5
    @test reciprocal(0) == 0.0
    @test reciprocal(2) == 1/2
    @test reciprocal_sum([1,2,5,0,0,4]) == (1/1 + 1/2 + 1/5 + 1/4)
end


@testset "missing utils" begin
    @test nonmissing(fill(missing, 5)) == 5
    @test nonmissing([1,2,3,missing, missing]) == 2
    @test nonmissing(x, "fca8") == 217
    @test nonmissings([1,2,missing, missing, 3], [3,missing, missing, 4, 5]) == [1,5]
end


end