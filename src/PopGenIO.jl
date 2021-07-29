module PopGenIO

# Write your package code here.
using CSV, DataFrames, NaturalSort, PooledArrays, Requires

include("Datasets.jl")
include("Delimited.jl")
include("Genepop.jl")
include("Read.jl")
include("Structure.jl")
include("PopData.jl")
include("Utils.jl")
include("VariantCall.jl")
@init @require  VariantCallFormat="28eba6e3-a997-4ad9-87c6-d933b8bca6c1" begin
    include("io/VariantCallLazy.jl")
end
@init @require VariantCallFormat="28eba6e3-a997-4ad9-87c6-d933b8bca6c1" begin
    @require GZip="92fee26a-97fe-5a0c-ad85-20a5f3185b63" include("io/VariantCallGzLazy.jl")
end

end
