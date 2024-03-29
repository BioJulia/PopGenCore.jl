"""
    isbinary(filepath::String)
Returns `true` if the `filepath` is a binary file. 
"""
function isbinary(filepath::String)
    !isvalid(String(Base.read(filepath, 1000)))
end

## Utilities for file reading/writing ##

#TODO change in docs
"""
    findploidy(genotypes::T where T<:AbstractVector)
Used internally in the `genepop` and `delimited` file parsers to scan the genotypes
of a sample and return the ploidy of the first non-missing locus.
"""
@inline function findploidy(genotypes::T) where T<:AbstractVector
    @inbounds Int8(length(first(skipmissing(genotypes))))
end


"""
    phase(loc::Union{String, Int}, type::DataType, digit::Int)
Takes a String of numbers or Integers and returns a typed locus appropriate for PopGen.jl as used in the
`genepop` and `csv` file parsers. Use `type` to specify output type (`Int8` or `Int16`),
and `digit` to specify the number of digits/characters used per allele in a locus.

## Examples
```
ph_locus = phase("128114", Int16, 3)
map(i -> phase(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```
"""
@inline function phase(loc::T, type::DataType, digit::Int) where T<:AbstractString
    loc == "-9" || iszero(parse(Int, loc)) && return missing
    phased = map(i -> parse(type, join(i)), Iterators.partition(loc, digit))
    NTuple{length(phased), type}(sort(phased))
end
phase(loc::Missing, type::DataType, digit::Int) = missing

@inline function phase(loc::T, type::DataType, digits::T) where T<:Integer
    loc == -9 || iszero(loc) && return missing
    out = type[]
    units = 10^digits
    d,r = divrem(loc, units)
    out = type[r]
    #@inbounds push!(out, r)
    @inbounds while d != 0
        d,r = divrem(d, units)
        @inbounds push!(out, r)
    end
    return NTuple{length(out), type}(sort(out))
end
precompile(phase, (Missing, DataType, Int64))
precompile(phase, (String, DataType, Int64))
precompile(phase, (String, DataType, Int64))
precompile(phase, (Int64, DataType, Int64))

"""
    unphase(geno::T; digits::Int = 3, ploidy::Int = 2, miss::Int = 0) where T <: Genotype
Takes a `Genotype` (e.g. `(131, 94)`) and returns a string of concatenated
alleles padded with *n* number of zeroes, where *n* is given by `digits = `.
`missing` values are returned as either a string of 'digits × ploidy' zeroes (`miss = 0`)
or `"-9"` (`miss = -9`). The `ploidy` flag is only relevant for unphasing `missing` genotypes
and not used otherwise.

### Example
```
unphase((1,2,3,4), digits = 3)
"001002003004"

unphase(missing, digits = 2, ploidy = 2, miss = -9)
"-9"

unphase(missing, digits = 2, ploidy = 2, miss = 0)
"0000"
```
"""
function unphase(geno::T; digits::Int = 3, ploidy::Int = 2, miss::Int = 0) where T <: Genotype
    join(map(i -> lpad(i, digits, "0"), geno))
end


function unphase(geno::Missing; digits::Int = 3, ploidy::Int, miss::Int = 0)
    if miss == 0
        return "0"^(digits * ploidy)
    else
        return "$miss"
    end
end
precompile(unphase, (NTuple{2, Int8},))
precompile(unphase, (NTuple{2, Int16},))
precompile(unphase, (Missing,))
