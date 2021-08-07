pfile = normpath(joinpath(@__DIR__, "..", "precompile", "precompile"))
genepop(pfile * ".gen", silent = true);
tst = structure(pfile * ".str", silent = true);