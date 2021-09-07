genepop(normpath(joinpath(@__DIR__, "..", "precompile", "precompile")) * ".gen", silent = true);
structure(normpath(joinpath(@__DIR__, "..", "precompile", "precompile")) * ".str", silent = true);
delimited(normpath(joinpath(@__DIR__, "..", "precompile", "precompile")) * ".csv", silent = true);
vcf(normpath(joinpath(@__DIR__, "..", "precompile", "precompile")) * ".vcf", silent = true) ;