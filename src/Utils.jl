export copy, drop_monomorphic, drop_monomorphic!, drop_multiallelic, drop_multiallelic!, loci, samples

#=
This file contains the helper functions necessary for file import/export
of various file formats.
=#


function Base.copy(data::PopData)
    PopData(copy(data.metadata), copy(data.genodata))
end


"""
    determine_marker(geno_parse::T, digits::Int) where T<:AbstractDataFrame
Return either `Int8` or `Int16` depending on largest allelic value in all genotypes
in the first 10 samples of an input DataFrame (or all the samples if less than 10 samples).
If the largest allele is 20 or greater, the marker will be considered a Microsatellite
and coded in `PopData` as `Int16`, and the opposite is true for SNPs. There's no
specific reason 100 was chosen other than it being a reasonable buffer for edge
cases since SNP data <= 4, and haplotyped data could be a bit higher. Even if the
microsatellite markers are coded incorrectly, there will be zero impact to performance,
and considering how few microsatellite markers are used in typical studies, the
effect on in-memory size should be negligible (as compared to SNPs).
"""
function determine_marker(geno_parse::T, digits::Int) where T<:AbstractDataFrame
    # get the total # columns
    total_col = size(geno_parse,2)
    # find the 25% cutoff
    if total_col > 11
        num_test_loc = (total_col - 2) ÷ 8
    else
        num_test_loc = total_col - 2
    end
    # remove everything else
    test_df = @view geno_parse[:, 3:(2 + num_test_loc)]

    # isolate the largest allele value
    max_allele = map(eachcol(test_df)) do i
        phase.(i, Int16, digits)  |>
        skipmissing |> Base.Iterators.flatten |> maximum
    end |> maximum

    if max_allele <= 100
        return Int8
    else
        return Int16
    end
end

"""
    find_ploidy(genotypes::T where T<:SubArray)
Used internally in the `genepop` and `delimited` file parsers to scan the genotypes
of a sample and return the ploidy of the first non-missing locus.
"""
@inline function find_ploidy(genotypes::T) where T<:GenoArray
    return Int8(length(first(skipmissing(genotypes))))
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
    sort!(phased)
    Tuple(phased)
end

phase(loc::Missing, type::DataType, digit::Int) = missing

@inline function phase(loc::T, type::DataType, digits::T) where T<:Integer
    loc == -9 || iszero(loc) && return missing
    out = type[]
    units = 10^digits
    d,r = divrem(loc, units)
    @inbounds push!(out, r)
    @inbounds while d != 0
        d,r = divrem(d, units)
        @inbounds push!(out, r)
    end
    return Tuple(sort(out))
end


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


function convert_coord(coordinate::String)
    lowercase(coordinate) == "missing" && return missing
    coord_strip = replace(uppercase(coordinate), r"[NSEW]" => "")
    split_coord = parse.(Float32, split(coord_strip, r"\s|:"))
    split_coord[2] /= 60.0
    if length(split_coord) == 3
        split_coord[3] /= 3600.0
    end
    conv = mapreduce(abs, +, split_coord)
    # N + E are positive | S + W are negative
    if split_coord[1] < 0 || occursin(r"[SW]", uppercase(coordinate))
        # negative
        return round(conv * -1, digits = 4)
    else
        # positive
        return round(conv, digits = 4)
    end
end

"""
    drop_monomorphic(data::PopData; silent::Bool = false)
Return a `PopData` object omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function drop_monomorphic(data::PopData; silent::Bool = false)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    monomorphs = string.(all_loci[[length(unique(i))==1 for i in prt]])
    if length(monomorphs) == 0
        return data
    elseif !silent
        if length(monomorphs) == 1
            println("Removing monomorphic locus: " * monomorphs[1])
            println()
        else
            println("Removing $(length(monomorphs)) monomorphic loci:" * "\n $monomorphs")
            println()
        end
    end
    exclude(data, locus = loci_to_rm)
end


"""
    drop_monomorphic!(data::PopData; silent::Bool = false)
Edit a `PopData` object in place by omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function drop_monomorphic!(data::PopData; silent::Bool = false)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    monomorphs = string.(all_loci[[length(unique(i))==1 for i in prt]])
    if length(monomorphs) == 0
        return data
    elseif !silent
        if length(monomorphs) == 1
            println("Removing monomorphic locus: " * monomorphs[1])
            println()
        else
            println("Removing $(length(monomorphs)) monomorphic loci:" * "\n $monomorphs")
            println()
        end
    end
    exclude!(data, locus = monomorphs)
end


"""
    drop_multiallelic(data::PopData)
Return a `PopData` object omitting loci that are not biallelic.
"""
function drop_multiallelic(data::PopData)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 0
        return data
    elseif length(nonbi) == 1
        @info "Removing 1 multiallelic locus"
        println()
    else
        @info "Removing $(length(nonbi)) multialleic loci"
        println()
    end
    exclude(data, locus = nonbi)
end


"""
    drop_multiallelic!(data::PopData)
Edit a `PopData` object in place, removing loci that are not biallelic.
"""
function drop_multiallelic!(data::PopData)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 0
        return data
    elseif length(nonbi) == 1
        @info "Removing 1 multiallelic locus"
        println()
    else
        @info "Removing $(length(nonbi)) multialleic loci"
        println()
    end
    exclude!(data, locus = nonbi)
end


"""
    generate_meta(data::DataFrame)
Given a genotype DataFrame formatted like `PopData.genodata`, generates a corresponding
`meta` DataFrame. In other words, it creates the `.metadata` part of `PopData` from the `.genodata` part.

**Example**
```
julia> cats = @nancycats ;

julia> cats_nometa = cats.genodata ;

julia> cats_meta = generate_meta(cats_nometa)
237×5 DataFrame
 Row │ name    population  ploidy  longitude  latitude 
     │ String  String      Int8    Float32?   Float32? 
─────┼─────────────────────────────────────────────────
   1 │ N215    1                2   missing   missing  
   2 │ N216    1                2   missing   missing  
   3 │ N217    1                2   missing   missing  
   4 │ N218    1                2   missing   missing  
   5 │ N219    1                2   missing   missing  
   6 │ N220    1                2   missing   missing  
   7 │ N221    1                2   missing   missing  
  ⋮  │   ⋮         ⋮         ⋮         ⋮         ⋮
 232 │ N295    17               2   missing   missing  
 233 │ N296    17               2   missing   missing  
 234 │ N297    17               2   missing   missing  
 235 │ N281    17               2   missing   missing  
 236 │ N289    17               2   missing   missing  
 237 │ N290    17               2   missing   missing  
                                       224 rows omitted
```
"""
function generate_meta(data::DataFrame)
    grp = groupby(data, :name)
    nms = map(z -> z.name, keys(grp))
    pops = [first(z.population) for z in grp]
    ploids = [find_ploidy(z.genotype) for z in grp]
    DataFrame(
        :name => nms,
        :population => pops,
        :ploidy => ploids,
        :longitude => Vector{Union{Missing, Float32}}(undef, (length(nms))),
        :latitude => Vector{Union{Missing, Float32}}(undef, (length(nms))),
        copycols = true
    )
end

"""
    Base.sort(x::NTuple{N,T}) where N where T <: Signed 
Sort the integers within a Tuple and return the sorted Tuple.
"""
function Base.sort(x::NTuple{N,T}) where N where T <: Signed 
    Tuple(sort(SVector(x)))
end

"""
    reciprocal(num::T) where T <: Signed
Returns the reciprocal (1/number) of a number. Will return `0` when
the number is `0` instead of returning `Inf`.
"""
function reciprocal(num::T) where T <: Real
    !iszero(num) ? 1.0/float(num) : 0.0
end


"""
    reciprocal_sum(x::AbstractVector{T}) where T<:Real
Return the sum of the reciprocal values of `x`, skipping the `Inf` values
resulting from divide-by-zero errors.
"""
function reciprocal_sum(x::AbstractVector{T}) where T<:Real
    mapreduce(reciprocal, +, x)
end

"""
    partitionarray(array::AbstractArray, steps::AbstractVector{<:Integer})
Like Base.Iterators.Partition, except you can apply arbitrary sizes to
partition the array by. The `steps` must add up to the total row length
of the array.

**Example**
```
julia> partitionmatrix(rand(20,5), [10,3,4,3]) .|> size
((10, 5), (3, 5), (4, 5), (3, 5))
```
"""
# solution brilliantly provided by @stevengj and @mcabbott on Slack and Discourse (https://discourse.julialang.org/t/is-there-a-simple-intuitive-way-to-partition-a-matrix-by-arbitrary-strides-like-i/55863)
function partitionarray(array::AbstractArray, steps::AbstractVector{<:Integer})
    v = axes(array,1)
    v == 1:sum(steps) || error("Steps provided do not sum to length of the first dimension")
    i = firstindex(v)
    tmp = (view(v, i:(i+=s)-1) for s in steps)
    [view(array,r,:) for r in tmp]
end

#=
"""
    loci(data::PopData)
Returns an array of strings of the loci names in a `PopData` object.
"""
loci(data::PopData) = data.genodata.locus.pool

"""
    samples(data::PopData)
View individual/sample names in a `PopData`
"""
samples(data::PopData) = @view data.metadata[!, :name]
=#