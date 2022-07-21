module  TestAlleleMatrices

using PopGenCore
using Test

cats = @nancycats ;

@testset "AlleleMatrices.jl" begin
    @testset "matrix" begin
        @test size(matrix(cats,"count", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(matrix(cats,"frequency", missings = "missing", scale = false, center = false)) == (237, 108)
        @test size(matrix(cats,"frequency", missings = "zero", scale = false, center = false)) == (237, 108)
        @test size(matrix(cats,"frequency", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(matrix(cats, scale = true)) == (237, 108)
        @test size(matrix(cats, center = true)) == (237, 108)
        @test matrix(cats, scale = true) != matrix(cats, center = true)
        a = matrix(cats, "count")
        b = matrix(cats, missings = "missing")
        c = matrix(cats, missings = "zero")
        d = matrix(cats)
        @test all(a .== b) == false
        @test all(a .== c) == false
        @test all(a .== d) == false
        @test all(b .=== c) == false
        @test all(b .=== d) == false
        @test all(c .== d) == false
    end
end

end # module