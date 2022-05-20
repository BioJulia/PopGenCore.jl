module  TestAlleleMatrices

using PopGenCore
using Test

cats = @nancycats ;

@testset "AlleleMatrices.jl" begin
    @testset "allelematrix" begin
        @test size(PopGen.allelematrix(cats, by = "count", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(PopGen.allelematrix(cats, by = "frequency", missings = "missing", scale = false, center = false)) == (237, 108)
        @test size(PopGen.allelematrix(cats, by = "frequency", missings = "zero", scale = false, center = false)) == (237, 108)
        @test size(PopGen.allelematrix(cats, by = "frequency", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(PopGen.allelematrix(cats, scale = true)) == (237, 108)
        @test size(PopGen.allelematrix(cats, center = true)) == (237, 108)
        @test PopGen.allelematrix(cats, scale = true) != PopGen.allelematrix(cats, center = true)
        a = PopGen.allelematrix(cats, by = "count")
        b = PopGen.allelematrix(cats, missings = "missing")
        c = PopGen.allelematrix(cats, missings = "zero")
        d = PopGen.allelematrix(cats)
        @test all(a .== b) == false
        @test all(a .== c) == false
        @test all(a .== d) == false
        @test all(b .=== c) == false
        @test all(b .=== d) == false
        @test all(c .== d) == false
    end
end

end # module