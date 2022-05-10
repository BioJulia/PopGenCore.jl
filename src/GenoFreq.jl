"""
    genocount_observed(locus::GenoArray)
Return a `Dict` of genotype counts of a single locus in a
`PopData` object.
"""
@inline function genocount_observed(locus::T) where T<:GenoArray
    # conditional testing if all genos are missing
    isallmissing(locus) && return missing
    countmap(skipmissing(locus))
end
precompile(genocount_observed, (Vector{Union{NTuple{2, Int8}}},))
precompile(genocount_observed, (Vector{Union{NTuple{2, Int16}}},))


#NOTICE the original method does not merge symmetrical genotypes
"""
    genocount_expected(locus::GenoArray)
Return a `Dict` of the expected genotype counts of a single locus in a
`PopData` object. Expected counts are calculated as the product of observed
allele frequencies multiplied by the number of non-missing genotypes.
"""
function genocount_expected(locus::T) where T<:GenoArray
    #count number of non-missing genotypes in the locus
    n = nonmissing(locus)
    ## get the observed allele frequencies
    allele_dict = allelefreq(locus)
    expected = Dict{nonmissingtype(eltype(T)), Float64}()
    @inbounds for (allele1, freq1) in allele_dict
        @inbounds for (allele2, freq2) in allele_dict
            geno = sort((allele1, allele2))
            #geno = (allele1, allele2)   # if not merging symmetrical genotypes
            expected[geno] = get!(expected, geno, 0.0) + (freq1 * freq2 * n)
        end
    end
    return expected
end
precompile(genocount_expected, (Vector{Union{NTuple{2, Int8}}},))
precompile(genocount_expected, (Vector{Union{NTuple{2, Int16}}},))


"""
    genofreq(locus::GenoArray)
Return a `Dict` of genotype frequencies of a single locus in a
`PopData` object.
"""
@inline function genofreq(locus::T) where T<:GenoArray
    # conditional testing if all genos are missing
    isallmissing(locus) && return missing
    proportionmap(locus |> skipmissing |> collect)
end
precompile(genofreq, (Vector{Union{NTuple{2, Int8}}},))
precompile(genofreq, (Vector{Union{NTuple{2, Int16}}},))


"""
    genofreq(data::PopData, locus::String; population::Bool = false)
Return a `Dict` of genotype frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of genotype frequencies
of that locus per population.
### Example
```
cats = @nancycats;
genofreq(cats, "fca8")
genofreq(cats, "fca8", population = true)
```
"""
function genofreq(data::PopData, locus::String; population::Bool=false)
    if !population
        tmp = data.genodata[data.genodata.locus .== locus, :genotype]
        genofreq(tmp)
    else
        tmp = data.genodata[data.genodata.locus .== locus, :]
        DataFrames.combine(groupby(tmp, :population), :genotype => genofreq => :freq)
    end
end
precompile(genofreq, (PopData, String))

"""
    genofreq_expected(locus::GenoArray)
Return a `Dict` of the expected genotype frequencies of a single locus in a
`PopData` object. Expected frequencies are calculated as the product of
observed allele frequencies.
"""
function genofreq_expected(locus::T) where T<:GenoArray
    ## get the observed allele frequencies
    allele_dict = allelefreq(locus)
    expected = Dict{nonmissingtype(eltype(T)), Float64}()
    @inbounds for (allele1, freq1) in allele_dict
        @inbounds for (allele2, freq2) in allele_dict
            geno = sort((allele1, allele2))
            #geno = (allele1, allele2)
            expected[geno] = get!(expected, geno, 0.0) + (freq1 * freq2)
        end
    end
    return expected
end
precompile(genofreq_expected, (Vector{NTuple{2, Int8}},))
precompile(genofreq_expected, (Vector{NTuple{2, Int16}},))



"""
    genofreq_expected(data::PopData, locus::String; population::Bool = false)
Return a `Dict` of expected genotype frequencies of a single locus in a
`PopData` object. Use `population = true` to return a table of expected genotype
frequencies of that locus per population.
### Example
```
cats = @nancycats;
genofreq_expected(cats, "fca8")
genofreq_expected(cats, "fca8", population = true)
```
"""
function genofreq_expected(data::PopData, locus::String; population::Bool=false)
    if !population
        tmp = data.genodata[data.genodata.locus .== locus, :genotype]
        genofreq_expected(tmp)
    else
        tmp = data.genodata[data.genodata.locus .== locus, :]
        DataFrames.combine(groupby(tmp, :population), :genotype => genofreq_expected => :freq)
    end
end
precompile(genofreq_expected, (PopData, String))