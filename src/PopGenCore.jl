module PopGenCore

using CSV, DataFrames, NaturalSort, PooledArrays, Requires, StaticArrays
using StatsBase: countmap

include("PopData.jl")
# Utilities
include("Utils/GeneralUtils.jl")
include("Utils/GenotypeUtils.jl")
include("Utils/ioUtils.jl")
include("Utils/MathUtils.jl")
include("Utils/MissingUtils.jl")
##
include("Conditionals.jl")
include("Permutations.jl")
include("Iterators.jl")
include("Manipulate.jl")
include("Delimited.jl")
include("Genepop.jl")
include("Structure.jl")
include("VariantCall.jl")
include("Read.jl")
include("Datasets.jl")
include("Permutations.jl")
@init @require  VariantCallFormat="28eba6e3-a997-4ad9-87c6-d933b8bca6c1" begin
    include("VariantCallLazy.jl")
end
@init @require VariantCallFormat="28eba6e3-a997-4ad9-87c6-d933b8bca6c1" begin
    @require GZip="92fee26a-97fe-5a0c-ad85-20a5f3185b63" include("VariantCallGzLazy.jl")
end

# precompile some file IO
include("precompile/precompile.jl") ;

end
