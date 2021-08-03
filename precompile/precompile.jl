pfile = normpath(joinpath(@__DIR__, "..", "precompile", "precompile"))
genepop(pfile * ".gen", silent = true);
structure(pfile * ".str", silent = true);