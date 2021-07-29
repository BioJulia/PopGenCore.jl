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
include("VariantCallGzLazy.jl")
include("VariantCallLazy.jl")

end
