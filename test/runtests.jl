fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
anyerrors = 0

using PopGenCore
using Test

all_tests = [
    "popdata.jl",
    "types.jl",
    "allelegenofreq.jl",
    "io.jl",
    "manipulate.jl",
    "utils.jl",
    "conditionals.jl",
    "iterators.jl"
]

println("Running tests:")

@testset "All Tests" begin
    for a_test in all_tests
        try
            include(a_test)
            println("\t\033[1m\033[32mPASSED\033[0m: $(a_test)")
        catch e
            global anyerrors += 1
            println("\t\033[1m\033[31mFAILED\033[0m: $(a_test)")
            if fatalerrors
                rethrow(e)
            elseif !quiet
                showerror(stdout, e, backtrace())
                println()
            end
        end
    end
end

anyerrors > 0 && throw("$anyerrors files of tests failed :(")
