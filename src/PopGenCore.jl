module PopGenCore

using CSV, DataFrames, NaturalSort, PooledArrays, Requires, StaticArrays
using StatsBase: countmap
using Random: shuffle, shuffle!

include("PopData.jl")
## Utilities
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
## File IO
include("io/Delimited.jl")
include("io/Genepop.jl")
include("io/Structure.jl")
include("io/VariantCall.jl")
include("io/Read.jl")
##
include("Datasets.jl")
@init @require  VariantCallFormat="28eba6e3-a997-4ad9-87c6-d933b8bca6c1" begin
    include("io/VariantCallLazy.jl")
end
@init @require VariantCallFormat="28eba6e3-a997-4ad9-87c6-d933b8bca6c1" begin
    @require GZip="92fee26a-97fe-5a0c-ad85-20a5f3185b63" include("io/VariantCallGzLazy.jl")
end

# precompile some file IO
include("precompile/precompile.jl") ;

end
