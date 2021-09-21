module  TestIterators

using PopGenCore
using Test

x = @nancycats

@testset "Iterators.jl" begin
    @testset "partition array" begin
        @test partitionarray(rand(20,5), [10,3,4,3]) .|> size == [(10, 5), (3, 5), (4, 5), (3, 5)]
        @test partitionarray(rand(20), [10,3,4,3]) .|> size == [(10, 1), (3, 1), (4, 1), (3, 1)]
    end

    @testset "pairiwse pairs" begin
        @test length(collect(pairwise_pairs(1:3))) == 3
        @test collect(pairwise_pairs(1:3)) == [(1, 2), (1, 3), (2, 3)]
        @test length(collect(pairwise_pairs(1:10))) == 45
        @test length(collect(pairwise_pairs(string.(1:3)))) == 3
        @test collect(pairwise_pairs(string.(1:3))) == [("1", "2"), ("1", "3"), ("2", "3")]
        @test length(collect(pairwise_pairs(string.(1:10)))) == 45
    end

    #TODO moved to PopGenSims.jl
    @testset "simulated pairs" begin
        @test sim_pairs(string.(1:4)) == [("1","2"), ("3","4")]
        @test_throws ArgumentError sim_pairs(string.(1:3))
    end

    @testset "skip iterators" begin
        @test collect(skipinf([1,2,3, Inf, -Inf, 4])) == [1,2,3,4]
        @test collect(skipnan([1,2,3, NaN,4])) == [1,2,3,4]
        @test collect(skipinfnan([1,2,3, Inf, -Inf, NaN, 4])) == [1,2,3,4]
    end
end

end